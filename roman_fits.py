"""FITS + DS9 output layer for roman-mast-tools.

Builds on the streaming primitives in `roman_mast`. All entry points here
take a streamed `{sca: AsdfFile}` dict (i.e. what `DataResults.stream()`
returns) plus the originating `Exposure`, so streaming stays decoupled from
output — you can stream once and dispatch to disk *and* DS9, or reuse the
same in-memory data for something else entirely.

Public surface
--------------
    to_fits_files(af_dict, exposure, out_dir=None, compress=False, ...)
        Write one FITS file per SCA, with SIP-approximated WCS headers.

    to_ds9(af_dict, exposure, dq_overlay=True, catalog_paths=None, ...)
        Pipe every SCA into a single DS9 WCS mosaic frame via XPA, optionally
        with a DQ mask overlay and/or L4 catalog source overlays. No local
        files touched (except the temp parquet cache for catalogs).

    download_catalogs(exposure, missions, out_dir=None, ...)
        Fetch the per-SCA L4 catalog (`cat_sca`) parquet files for an
        exposure. Returns {sca: local_path or None}. Silently returns None
        for SCAs that don't have a catalog on MAST.

Both are also exposed as convenience methods `DataResults.to_fits()` and
`DataResults.to_ds9()` in `roman_mast`, which do the stream + output in
one call.

Notes on the SIP WCS
--------------------
Roman ships a gwcs object per SCA. DS9 / astropy consume FITS headers only,
so we approximate gwcs → FITS SIP via `wcs.to_fits_sip(degree=sip_degree)`.
Default degree=4 is accurate to ~0.1 pixel across a Roman SCA, which is
plenty for DS9 mosaicking and typical analysis.
"""

from __future__ import annotations

import io
import os
from typing import Any, Optional

import numpy as np
from astropy.io import fits
import roman_datamodels as rdm

from roman_mast import Exposure, _log, close_streams


# Which array attributes we pre-fetch off the datamodel. Kept to the two
# arrays the current sinks actually consume:
#   - `data` → FITS ImageHDU pixels, DS9 science mosaic, CSV data stats
#   - `dq`   → DS9 mask overlay
# Add `err` / `var_poisson` / `var_rnoise` / `var_flat` here if a future
# sink wants uncertainty maps — each one is another ~17-67 MB per SCA
# (×18 SCAs), so we leave them out until something needs them.
_MATERIALIZE_ATTRS = ('data', 'dq')


def _materialize_dm(dm) -> None:
    """Force-load a datamodel's big array attributes into in-memory ndarrays.

    Roman AsdfFiles stream lazily via fsspec + a MAST-issued pre-signed S3
    URL that expires after 60 s. If we return the dm and let a downstream
    sink read `dm.data` (or worse, `dm.data` *then* `dm.dq`, which live at
    different block offsets and re-fetch from S3) minutes later, that second
    fetch hits HTTP 403. Reading every array we care about here — while the
    URL is still fresh — makes the dm self-contained; sinks never re-hit S3.

    Writes materialized ndarrays back onto the dm in place when the model
    accepts it, so `dm.data` etc. keep working transparently.
    """
    for attr in _MATERIALIZE_ATTRS:
        try:
            arr = getattr(dm, attr, None)
        except (AttributeError, KeyError):
            continue
        if arr is None:
            continue
        try:
            materialized = np.asarray(arr[...])
        except Exception:
            continue  # Not an array-like attribute after all — skip.
        try:
            setattr(dm, attr, materialized)
        except Exception:
            # Datamodel refused the assignment (some fields are validated).
            # The read still populated fsspec's block cache, so subsequent
            # reads within the URL lifetime will hit the cache. This is
            # a best-effort belt-and-suspenders — the main win is the read.
            pass


def _stream_one_sca(sca: int, filename: str, missions):
    """Worker: open one SCA, wrap it, materialize its arrays. Thread-safe.

    Returns (sca, dm) on success, (sca, None) on failure. Failures are
    logged here so the pool doesn't swallow them silently.
    """
    try:
        af = missions.read_product(filename)
        dm = rdm.open(af)
        _materialize_dm(dm)
        return sca, dm
    except Exception as e:
        _log(f"ERROR streaming SCA {sca:02d} ({filename}): "
             f"{type(e).__name__}: {e}")
        return sca, None


def stream_materialized(exposure: Exposure, missions, *, scas=None,
                        show_progress=True, max_workers: int = 8) -> dict:
    """Stream every SCA of `exposure` in parallel and pre-load arrays inline.

    Like `roman_mast.stream_exposure`, but (a) interleaves the S3 open with
    an immediate array read per SCA — so no consumer downstream ever needs
    the (60 s) pre-signed URL again — and (b) fans the per-SCA reads out
    over a thread pool. Threads are the right shape here: each `read_product`
    call is a blocking-network-I/O read that releases the GIL, and MAST
    mints an independent pre-signed URL per request so there's nothing
    shared for the threads to contend on.

    Parameters
    ----------
    max_workers : int
        Concurrent SCA fetches. Default 8 — empirically saturates the
        per-connection bandwidth on the /home/mrizzo node without tripping
        S3 request throttling. Set to 1 to force sequential fetching if a
        machine is bandwidth-shared or you want deterministic ordering in
        logs.

    Returns
    -------
    dict[int, DataModel]
        Mapping SCA integer → open roman_datamodels DataModel with its
        big arrays already resident in memory. A None value marks an SCA
        we failed to stream.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    try:
        from tqdm.auto import tqdm
    except ImportError:  # pragma: no cover
        def tqdm(iterable, **_):
            return iterable

    if scas is not None:
        wanted = set(int(s) for s in scas)
        pairs = [(s, f) for s, f in zip(exposure.scas, exposure.filenames)
                 if s in wanted]
        missing = wanted - {s for s, _ in pairs}
        if missing:
            _log(f"WARNING: exposure has no filenames for SCAs {sorted(missing)}")
    else:
        pairs = list(zip(exposure.scas, exposure.filenames))

    workers = max(1, min(max_workers, len(pairs)))
    _log(f"Streaming exposure visit_id={exposure.visit_id} "
         f"exp={exposure.exposure} ({len(pairs)} SCA files, "
         f"{workers} workers, materializing inline)")

    dm_dict: dict = {}

    if workers == 1:
        # Preserve the sequential codepath verbatim — one less thing to
        # debug when a machine misbehaves under threading.
        iterator = tqdm(pairs, desc="Streaming SCAs", disable=not show_progress)
        for sca, filename in iterator:
            _, dm = _stream_one_sca(sca, filename, missions)
            dm_dict[sca] = dm
        return dm_dict

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(_stream_one_sca, sca, filename, missions)
            for sca, filename in pairs
        ]
        iterator = tqdm(
            as_completed(futures), total=len(futures),
            desc="Streaming SCAs", disable=not show_progress,
        )
        for fut in iterator:
            sca, dm = fut.result()
            dm_dict[sca] = dm

    # Return in SCA order so downstream logs/CSV rows come out ascending.
    return {sca: dm_dict[sca] for sca in sorted(dm_dict)}


# ---------------------------------------------------------------------------
# Optional pyds9 — only needed for to_ds9()
# ---------------------------------------------------------------------------

try:
    import pyds9
except ImportError:  # pragma: no cover
    pyds9 = None


# ---------------------------------------------------------------------------
# Shared header builder
# ---------------------------------------------------------------------------

# meta path (dot-separated) → (FITS keyword, comment)
_META_KEYWORDS = [
    ('exposure.start_time',              'DATE-BEG', 'Exposure start (UTC)'),
    ('exposure.end_time',                'DATE-END', 'Exposure end (UTC)'),
    ('exposure.effective_exposure_time', 'EXPTIME',  'Effective exposure time [s]'),
    ('instrument.detector',              'DETECTOR', 'Detector name'),
    ('instrument.optical_element',       'FILTER',   'Optical element / filter'),
]


def _build_sca_header(dm, exposure: Exposure, scanum: int, sip_degree: int) -> fits.Header:
    """Return a FITS header for one SCA — SIP WCS + standard Roman keywords.

    Header is written such that both the plain ImageHDU (disk mode) and the
    mosaic-tile ImageHDU (DS9 mode) can use it verbatim; low-level keywords
    (SIMPLE/BITPIX/NAXIS*) are set by astropy from the data at HDU build time.
    """
    wcs = dm.meta.wcs
    hdr = wcs.to_fits_sip(bounding_box=wcs.bounding_box, degree=sip_degree)

    hdr['EXTNAME'] = f'SCA{scanum:02d}'
    hdr['VISITID'] = (exposure.visit_id, 'Roman visit ID')
    hdr['EXPNUM']  = (exposure.exposure, 'Exposure number within visit')
    hdr['SCANUM']  = (scanum, 'SCA number (1-18)')

    for src, dest, comment in _META_KEYWORDS:
        try:
            val = dm.meta
            for attr in src.split('.'):
                val = getattr(val, attr)
            hdr[dest] = (str(val), comment)
        except AttributeError:
            pass

    return hdr


def _open_dm(af):
    """Convert an AsdfFile → roman_datamodels DataModel. Passthrough if already one."""
    if hasattr(af, 'meta') and hasattr(af, 'data'):
        return af
    return rdm.open(af)


# ---------------------------------------------------------------------------
# to_fits_files — one SCA per FITS file on disk
# ---------------------------------------------------------------------------

def to_fits_files(
    af_dict: dict,
    exposure: Exposure,
    *,
    out_dir: Optional[str] = None,
    compress: bool = False,
    sip_degree: int = 4,
    overwrite: bool = True,
) -> str:
    """Write each streamed SCA to its own FITS file, one per SCA.

    Parameters
    ----------
    af_dict : dict[int, AsdfFile]
        The `{sca: AsdfFile}` mapping returned by `DataResults.stream()` /
        `stream_exposure()`.
    exposure : Exposure
        The exposure these SCAs come from (used for header keywords + the
        default output directory name).
    out_dir : str, optional
        Where to write. Defaults to ``wfi_fits_{visit_id}_exp{NN}/`` in cwd.
    compress : bool
        If True, write RICE_1 tile-compressed ``.fits.fz`` files (~3–4× smaller
        for background-limited data). Default False.
    sip_degree : int
        SIP polynomial degree for the gwcs → FITS approximation. Default 4.
    overwrite : bool
        If True, overwrite existing files. Default True.

    Returns
    -------
    str
        The output directory path.
    """
    if out_dir is None:
        out_dir = (f'wfi_fits_{exposure.visit_id}_exp{exposure.exposure:02d}')
    os.makedirs(out_dir, exist_ok=True)
    ext = 'fits.fz' if compress else 'fits'

    _log(f"Writing FITS files → {out_dir}/  (compress={compress}, sip_degree={sip_degree})")

    written = 0
    for scanum in sorted(af_dict):
        af = af_dict[scanum]
        if af is None:
            _log(f"  SCA {scanum:02d}: no data, skipping")
            continue

        dm = _open_dm(af)
        data = np.asarray(dm.data[...]).astype(np.float32)
        hdr = _build_sca_header(dm, exposure, scanum, sip_degree)

        if compress:
            comp = fits.CompImageHDU(
                data=data, header=hdr,
                compression_type='RICE_1', tile_shape=(256, 256),
            )
            hdul = fits.HDUList([fits.PrimaryHDU(), comp])
        else:
            hdul = fits.HDUList([fits.PrimaryHDU(data=data, header=hdr)])

        out_path = os.path.join(out_dir, f'sca_{scanum:02d}.{ext}')
        hdul.writeto(out_path, overwrite=overwrite)
        size_mb = os.path.getsize(out_path) / 1e6
        _log(f"  SCA {scanum:02d}: wrote {out_path}  "
             f"({data.shape[0]}x{data.shape[1]}, {size_mb:.0f} MB)")
        written += 1

    _log(f"Wrote {written} FITS file(s) in {out_dir}/")
    _log(f"  Open in DS9:  ds9 -mosaic {out_dir}/sca_*.{ext}")
    return out_dir


# ---------------------------------------------------------------------------
# L4 catalog (cat_sca) download + DS9 region synthesis
# ---------------------------------------------------------------------------

# `read_product` only handles .fits / .asdf; parquet has to go through
# `download_file`, which writes to disk. We cache into a per-run temp dir
# and clean up after the DS9 handoff.
def _catalog_filename(exposure: Exposure, scanum: int) -> str:
    """Build the expected `_cat.parquet` filename for one SCA of an exposure.

    Roman filenames follow the fixed pattern
    ``r{visit_id}_{exp:04d}_wfi{sca:02d}_{filter}_cat.parquet``. `filter`
    is the exposure's optical_element lowercased.
    """
    filt = (exposure.optical_element or 'unk').lower()
    return (f"r{exposure.visit_id}_{exposure.exposure:04d}_"
            f"wfi{scanum:02d}_{filt}_cat.parquet")


def download_catalogs(
    exposure: Exposure,
    missions,
    *,
    scas=None,
    out_dir: Optional[str] = None,
    max_workers: int = 8,
) -> dict:
    """Download L4 per-SCA catalogs for `exposure`. Returns {sca: path or None}.

    Not every SCA has a catalog on MAST — the L4 pipeline may skip an SCA
    if the calibrated image failed or lacked detections. `download_file`
    returns ``status='ERROR'`` with a 404 message in that case; we swallow
    those and return None for that SCA so callers can attempt the whole
    set unconditionally.

    Parameters
    ----------
    out_dir : str, optional
        Where to drop the parquet files. Default: a fresh tempdir (caller
        can shutil.rmtree it — the paths returned point inside it).
    max_workers : int
        Concurrent downloads. Matches the streaming pool; each catalog
        is a small (~50-200 kB) HTTP GET so parallelism helps latency
        but hits diminishing returns fast.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import tempfile

    if out_dir is None:
        out_dir = tempfile.mkdtemp(prefix='roman_cat_')
    os.makedirs(out_dir, exist_ok=True)

    if scas is None:
        wanted = list(exposure.scas)
    else:
        wanted = [s for s in exposure.scas if s in set(int(x) for x in scas)]

    def _fetch(sca):
        fname = _catalog_filename(exposure, sca)
        # download_file writes {out_dir}/{fname}. It returns (status, msg, url);
        # status is 'COMPLETE' on success, 'ERROR' on 404 / anything else.
        try:
            status, msg, _ = missions.download_file(
                fname, local_path=out_dir, cache=True, verbose=False,
            )
        except Exception as e:
            _log(f"  SCA {sca:02d}: catalog download crashed: "
                 f"{type(e).__name__}: {e}")
            return sca, None
        if status != 'COMPLETE':
            return sca, None
        return sca, os.path.join(out_dir, fname)

    workers = max(1, min(max_workers, len(wanted)))
    paths: dict = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_fetch, s) for s in wanted]
        for fut in as_completed(futures):
            sca, path = fut.result()
            paths[sca] = path

    have = sum(1 for p in paths.values() if p is not None)
    _log(f"L4 catalogs: {have}/{len(wanted)} SCAs have a cat_sca on MAST")
    return {s: paths[s] for s in sorted(paths)}


# Columns we pull off each catalog parquet to build the region file. Kept
# tight — pyarrow reads only these, so parquets stay cheap to load. Adding
# a column here also automatically makes it available for label formatting.
_REGION_COLUMNS = (
    'label',              # per-SCA integer source ID (unique inside one SCA)
    'ra', 'dec',          # sky coords for the region circle
    'warning_flags',      # 0 = clean; drop non-zero unless --catalog-include-flagged
    'kron_abmag',         # total mag — the number a user usually wants to see
    'kron_abmag_err',     # its uncertainty
    'is_extended',        # star (False) vs. galaxy-like (True); drives color
)


def _format_label(mode: str, sca: int, label, mag, mag_err) -> str:
    """Produce the DS9 `text={...}` payload for one source.

    `mode` mirrors --catalog-label:
        'none' → no text (returns '')
        'id'   → 'SCA01.42'
        'mag'  → '19.34'
        'full' → 'SCA01.42  m=19.34±0.02'
    Missing / NaN values render as '?' rather than 'nan' so a glance at DS9
    tells you what was in the catalog vs. what was blanked.
    """
    if mode == 'none':
        return ''
    sid = f"SCA{sca:02d}.{int(label)}" if label is not None else f"SCA{sca:02d}"
    if mode == 'id':
        return sid
    m = f"{mag:.2f}" if mag is not None and np.isfinite(mag) else '?'
    if mode == 'mag':
        return m
    e = f"{mag_err:.2f}" if mag_err is not None and np.isfinite(mag_err) else '?'
    return f"{sid}  m={m}±{e}"


# One region file (fk5 coordinates) covers the whole mosaic — DS9 evaluates
# each region against the currently displayed WCS, so we don't need per-SCA
# region files.
def _build_region_text(
    catalog_paths: dict,
    *,
    radius_arcsec: float = 0.4,
    color: str = 'green',
    extended_color: Optional[str] = 'yellow',
    width: int = 1,
    max_per_sca: Optional[int] = None,
    include_flagged: bool = False,
    label_mode: str = 'full',
) -> str:
    """Assemble a DS9 region string from a set of per-SCA catalog parquets.

    One `circle` per detected source, in fk5 sky coordinates. The mosaic
    WCS handles the projection, so a single region block works across all
    SCA tiles.

    Each region carries per-source metadata as DS9 comment properties so
    clicking a circle in DS9 shows the source ID / magnitude:
        text={SCA01.42  m=19.34±0.02}   — label shown next to the circle
        tag={SCA01} tag={extended}      — grouping tags (toggle per-SCA / per-class)

    Parameters
    ----------
    color : str
        Default region color; also used for stars (is_extended=False).
    extended_color : str, optional
        Color for `is_extended=True` sources. Set None to disable and use
        `color` for everything.
    max_per_sca : int, optional
        Cap on regions per SCA — keeps DS9 responsive for very rich fields.
        None (default) draws every source.
    include_flagged : bool
        By default we drop sources with a non-zero warning_flags value
        (saturation, edge, contamination, ...). Flip to True to show
        everything.
    label_mode : {'none', 'id', 'mag', 'full'}
        What to write into each region's DS9 label. 'none' = no text
        (fastest to render), 'id' = 'SCA01.42', 'mag' = '19.34',
        'full' = 'SCA01.42  m=19.34±0.02'. Default 'id'.
    """
    import pyarrow.parquet as pq

    lines = ['# DS9 region file — roman-mast-tools L4 sources',
             'global color=%s width=%d font="helvetica 10 normal"' % (color, width),
             'fk5']

    total = 0
    for scanum in sorted(catalog_paths):
        path = catalog_paths[scanum]
        if path is None:
            continue
        try:
            tbl = pq.read_table(path, columns=list(_REGION_COLUMNS))
        except Exception as e:
            _log(f"  SCA {scanum:02d}: could not read catalog {path}: "
                 f"{type(e).__name__}: {e}")
            continue

        # to_numpy() on nullable columns returns object arrays; that's fine
        # here — we only index element-wise, no vectorized arithmetic.
        labels    = tbl.column('label').to_numpy()
        ras       = tbl.column('ra').to_numpy()
        decs      = tbl.column('dec').to_numpy()
        flags     = tbl.column('warning_flags').to_numpy()
        mags      = tbl.column('kron_abmag').to_numpy()
        mag_errs  = tbl.column('kron_abmag_err').to_numpy()
        extended  = tbl.column('is_extended').to_numpy()

        # Drop rows with masked / non-finite sky coords — the L4 pipeline
        # can emit sources whose centroid failed to converge.
        good = np.isfinite(ras.astype(float)) & np.isfinite(decs.astype(float))
        if not include_flagged:
            good &= (flags == 0)

        idx = np.where(good)[0]
        if max_per_sca is not None and len(idx) > max_per_sca:
            idx = idx[:max_per_sca]

        for i in idx:
            props = [f'tag={{SCA{scanum:02d}}}']

            # extended → different color + a class tag so users can toggle
            # 'show galaxies only' via DS9's region-group filter.
            is_ext = bool(extended[i]) if extended[i] is not None else False
            if is_ext:
                props.append('tag={extended}')
                if extended_color:
                    props.append(f'color={extended_color}')
            else:
                props.append('tag={point}')

            text = _format_label(label_mode, scanum, labels[i], mags[i], mag_errs[i])
            if text:
                # DS9 text is single-braced; escape any stray '}' just in case.
                safe = text.replace('}', ')')
                props.append(f'text={{{safe}}}')

            lines.append(
                f'circle({ras[i]:.7f},{decs[i]:.7f},{radius_arcsec:.2f}") # '
                + ' '.join(props)
            )
        total += len(idx)

    _log(f"L4 catalog overlay: built {total} source regions "
         f"(label_mode={label_mode})")
    return '\n'.join(lines) + '\n'


# ---------------------------------------------------------------------------
# to_ds9 — pipe an in-memory MEF via XPA
# ---------------------------------------------------------------------------

def to_ds9(
    af_dict: dict,
    exposure: Exposure,
    *,
    sip_degree: int = 4,
    dq_overlay: bool = True,
    catalog_paths: Optional[dict] = None,
    catalog_radius_arcsec: float = 0.4,
    catalog_color: str = 'green',
    catalog_extended_color: Optional[str] = 'yellow',
    catalog_include_flagged: bool = False,
    catalog_label_mode: str = 'full',
    catalog_show_labels: bool = False,
    ds9_target: Optional[str] = None,
):
    """Stream every SCA into a single DS9 WCS mosaic frame via XPA.

    Nothing hits the local disk. Two in-memory MEFs are built (one for the
    science data, one for a per-pixel bad-pixel mask) and shipped to DS9 with
    ``fits mosaicimage wcs`` — the same command DS9 uses for on-disk mosaics,
    just fed from a byte buffer.

    Parameters
    ----------
    af_dict, exposure
        As in `to_fits_files`.
    sip_degree : int
        SIP polynomial degree for the gwcs → FITS approximation. Default 4.
    dq_overlay : bool
        If True (default), send DQ (dq != 0) as a mask overlay on top of the
        data frame. Silently skipped for products without a DQ array.
    catalog_paths : dict, optional
        ``{sca: path_to_cat_sca.parquet or None}`` — typically what
        `download_catalogs()` returns. When provided, DS9 draws one circle
        per source in fk5 coordinates on top of the mosaic. None values
        (SCAs with no catalog on MAST) are silently skipped.
    catalog_radius_arcsec : float
        Circle radius for source overlays. Default 0.4" — big enough to
        spot in a wide zoom, small enough not to obscure the source itself.
    catalog_color : str
        DS9 color name for point-like sources (default 'green').
    catalog_extended_color : str, optional
        Color for sources flagged is_extended=True. Default 'yellow' — gives
        a stars vs. galaxies split at a glance. None disables the split and
        uses `catalog_color` for everything.
    catalog_include_flagged : bool
        If False (default), sources with non-zero `warning_flags` are dropped.
    catalog_label_mode : {'none', 'id', 'mag', 'full'}
        DS9 label attached to each region. 'full' → 'SCA01.42  m=19.34±0.02'
        (default), 'id' → 'SCA01.42', 'mag' → '19.34', 'none' → no label
        (fastest for very rich fields).
    catalog_show_labels : bool
        If False (default), send `regions showtext no` to DS9 after loading —
        labels stay in the region metadata (visible when you click a source,
        and preserved in the saved .reg file) but don't paint on the image.
        Flip to True for always-on labels.
    ds9_target : str, optional
        Name of an existing DS9 XPA target. If None, connects to the default
        DS9 instance (starting one if pyds9 does that on your setup).

    Returns
    -------
    pyds9.DS9
        The DS9 instance, with the exposure loaded and zoomed-to-fit.
    """
    if pyds9 is None:
        raise ImportError(
            "pyds9 is required for to_ds9. Install with `pip install pyds9` "
            "(and make sure DS9 itself is running: `ds9 &`)."
        )

    try:
        d = pyds9.DS9(target=ds9_target) if ds9_target else pyds9.DS9()
    except Exception as e:
        raise RuntimeError(
            f"Failed to connect to DS9: {e}\n"
            f"Make sure DS9 is running: ds9 &"
        )

    _log(f"Building in-memory MEFs for DS9 (visit={exposure.visit_id} "
         f"exp={exposure.exposure}, sip_degree={sip_degree})")

    data_hdulist = fits.HDUList([fits.PrimaryHDU()])
    dq_hdulist   = fits.HDUList([fits.PrimaryHDU()])

    n_data = 0
    n_dq   = 0

    for scanum in sorted(af_dict):
        af = af_dict[scanum]
        if af is None:
            _log(f"  SCA {scanum:02d}: no data, skipping")
            continue

        dm = _open_dm(af)
        data = np.asarray(dm.data[...]).astype(np.float32)
        hdr  = _build_sca_header(dm, exposure, scanum, sip_degree)

        data_hdulist.append(
            fits.ImageHDU(data=data, header=hdr.copy(), name=f'SCA{scanum:02d}')
        )
        n_data += 1

        if dq_overlay:
            try:
                dq = np.asarray(dm.dq[...])
                dq_mask = (dq != 0).astype(np.uint8)
                dq_hdr = hdr.copy()
                dq_hdr['BUNIT']   = 'flag'
                dq_hdr['CONTENT'] = ('DQ_MASK', 'Non-zero = bad pixel (any DQ bit set)')
                dq_hdulist.append(
                    fits.ImageHDU(data=dq_mask, header=dq_hdr, name=f'DQ{scanum:02d}')
                )
                n_dq += 1
            except AttributeError:
                pass  # L1 uncal has no DQ layer yet — fine.

        _log(f"  SCA {scanum:02d}: added to mosaic ({data.shape[0]}x{data.shape[1]})")

    if n_data == 0:
        raise RuntimeError("No SCA data streamed; nothing to load into DS9.")

    # --- Science mosaic ---------------------------------------------------
    data_buf = io.BytesIO()
    data_hdulist.writeto(data_buf)
    data_bytes = data_buf.getvalue()

    _log(f"Piping {len(data_bytes)/1e6:.0f} MB data MEF into DS9 as WCS mosaic")
    d.set('frame delete all')
    d.set('frame new')
    d.set('fits mosaicimage wcs', data_bytes)

    # --- DQ overlay -------------------------------------------------------
    if dq_overlay and n_dq > 0:
        dq_buf = io.BytesIO()
        dq_hdulist.writeto(dq_buf)
        dq_bytes = dq_buf.getvalue()

        _log(f"Piping {len(dq_bytes)/1e6:.1f} MB DQ mask MEF into DS9 as overlay")
        d.set('mask clear')
        d.set('mask color red')
        d.set('mask transparency 50')
        d.set('mask mark nonzero')
        d.set('fits mask mosaicimage wcs', dq_bytes)

    # --- L4 catalog source overlay ---------------------------------------
    # DS9 evaluates fk5 regions against the mosaic's WCS, so one region
    # block spans every SCA. `regions` via XPA takes the region text on
    # stdin — same as loading a .reg file.
    n_cat = 0
    if catalog_paths:
        region_text = _build_region_text(
            catalog_paths,
            radius_arcsec=catalog_radius_arcsec,
            color=catalog_color,
            extended_color=catalog_extended_color,
            include_flagged=catalog_include_flagged,
            label_mode=catalog_label_mode,
        )
        # Count SCAs with any usable catalog — for the log line below.
        n_cat = sum(1 for p in catalog_paths.values() if p is not None)
        try:
            # 'regions -format ds9' reads a full region file; we send it
            # by pipe (the second positional arg to pyds9.DS9.set).
            d.set('regions -format ds9', region_text)
            # Hide the text=… labels by default so 900+ mag readouts don't
            # obscure the image. The labels stay in the region metadata,
            # so clicking a source pops them in DS9's region-info dialog.
            if not catalog_show_labels:
                d.set('regions showtext no')
        except Exception as e:
            _log(f"WARNING: could not load catalog regions into DS9: "
                 f"{type(e).__name__}: {e}")
            n_cat = 0

    d.set('zoom to fit')

    extras = []
    if n_dq:
        extras.append(f"DQ overlay ({n_dq} SCAs)")
    if n_cat:
        extras.append(f"L4 catalog ({n_cat} SCAs)")
    extras_str = f" + {' + '.join(extras)}" if extras else ""
    _log(f"Loaded {n_data} SCAs into DS9{extras_str}  "
         f"(visit {exposure.visit_id} / exposure {exposure.exposure})")

    return d


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _resolve_exposure_keys(res, spec):
    """Turn --exposures spec (int / 'a,b' / 'a-b' / 'all') into ordered indices.

    Returns 1-based indices into res.exposures. Raises IndexError on bad values.
    """
    from roman_mast import parse_int_spec

    if spec is None or str(spec).strip().lower() in ('', 'all', '*'):
        return list(range(1, res.n_exposures + 1))

    indices = parse_int_spec(spec)
    for i in indices:
        if i < 1 or i > res.n_exposures:
            raise IndexError(
                f"Exposure index {i} out of range 1..{res.n_exposures} "
                f"(query returned {res.n_exposures} exposure(s))"
            )
    return indices


def _cli():
    import argparse
    import warnings
    warnings.filterwarnings('ignore')

    from roman_mast import (
        add_list_data_args, list_data_from_args, print_summary,
    )

    p = argparse.ArgumentParser(
        description="Stream Roman WFI SCAs into FITS files or DS9 (WCS mosaic).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Selection works in two steps: the standard --program / --pass / --detector /
--visit-id filters find the exposure(s), then --exposures picks which of
those to output. --list first, then re-run without --list.

Every run drops one folder per exposure under --out-dir (default cwd):

    <out-dir>/v{visit_id}_exp{NN}/
        metadata_{visit_id}_exp{NN}.csv    (unless --no-metadata)
        catalog/                            (ds9 mode, unless --no-catalog)
            r..._wfi{NN}_..._cat.parquet    (one per SCA that has a catalog)
            catalog_{visit_id}_exp{NN}.reg  (DS9 regions from all catalogs)
        sca_{NN}.fits                       (fits mode)

Examples:
  # See what's available (no output yet)
  python roman_fits.py --program 114 --pass 57 --sca-only --list

  # Write one exposure to disk as 18 FITS files under /tmp/wfi/v..._exp01/
  python roman_fits.py --program 114 --pass 57 --sca-only \\
      --exposures 1 --to fits --out-dir /tmp/wfi

  # Same, RICE compressed (folder in cwd)
  python roman_fits.py --program 114 --pass 57 --sca-only \\
      --exposures 1 --to fits --compress

  # Multiple exposures — one v..._expNN/ folder each under --out-dir
  python roman_fits.py --program 114 --pass 57 --sca-only \\
      --exposures 1-4 --to fits --out-dir /tmp/wfi

  # Stream to DS9 (needs `ds9 &` running and pyds9 installed); catalog
  # parquets + .reg land in <cwd>/v..._exp01/catalog/
  python roman_fits.py --program 114 --pass 57 --sca-only \\
      --exposures 1 --to ds9

  # Only some SCAs
  python roman_fits.py --program 114 --pass 57 --sca-only \\
      --exposures 1 --scas 1-6 --to fits
""",
    )

    add_list_data_args(p)

    p.add_argument('--exposures', default='1',
                   help="Which exposure(s) to output (1-based index into the "
                        "listed exposures). '1', '1,3,5', '1-4', or 'all'. "
                        "Default: '1'.")
    p.add_argument('--scas', default=None,
                   help="Restrict to a subset of SCAs, e.g. '4' or '1-6' or "
                        "'1,3,5'. Default: every SCA the exposure has.")
    p.add_argument('--to', choices=['fits', 'ds9'], default='fits',
                   help="Output mode. 'fits' writes per-SCA FITS files; "
                        "'ds9' streams into a running DS9 as a WCS mosaic. "
                        "Default: fits.")
    p.add_argument('--out-dir', default=None,
                   help="Root directory for all per-exposure exports. Each "
                        "exposure gets a subfolder v{visit_id}_exp{NN}/ under "
                        "this root, containing the metadata CSV, the catalog/ "
                        "subfolder (ds9 mode), and — in fits mode — the per-SCA "
                        "FITS files. Default: cwd.")
    p.add_argument('--compress', action='store_true',
                   help='RICE_1 tile-compressed .fits.fz (fits mode only)')
    p.add_argument('--sip-degree', type=int, default=4,
                   help='SIP polynomial degree for gwcs → FITS (default 4)')
    p.add_argument('--no-dq-overlay', action='store_true',
                   help='Skip the DQ (bad-pixel) mask overlay in ds9 mode')
    p.add_argument('--no-catalog', action='store_true',
                   help='Skip the L4 per-SCA catalog (cat_sca) source overlay '
                        'in ds9 mode. By default, if a catalog exists on MAST '
                        "for each SCA it's downloaded and drawn as regions.")
    p.add_argument('--catalog-radius', type=float, default=0.4,
                   help='DS9 region radius (arcsec) for catalog sources '
                        '(default 0.4)')
    p.add_argument('--catalog-color', default='green',
                   help='DS9 color for point-like catalog sources '
                        '(is_extended=False). Default green.')
    p.add_argument('--catalog-extended-color', default='yellow',
                   help='DS9 color for extended catalog sources '
                        '(is_extended=True). Set to "same" to use '
                        '--catalog-color for everything. Default yellow.')
    p.add_argument('--catalog-label', choices=['none', 'id', 'mag', 'full'],
                   default='full',
                   help="What to store in each region's text= field. 'none' = "
                        "empty, 'id' = SCA{NN}.{label}, 'mag' = Kron AB mag, "
                        "'full' = 'SCA{NN}.{label}  m={mag}±{err}' (default). "
                        "By default DS9 renders these only when a region is "
                        "selected — use --catalog-show-labels for always-on.")
    p.add_argument('--catalog-show-labels', action='store_true',
                   help='Always render region labels in DS9. Default: labels '
                        'stay in region metadata but hidden until a source '
                        'is clicked (via "regions showtext no").')
    p.add_argument('--catalog-include-flagged', action='store_true',
                   help='Include sources with non-zero warning_flags. Default '
                        'drops them (saturation / edge / contamination).')
    p.add_argument('--catalog-dir', default=None,
                   help='Override the per-exposure catalog folder location. '
                        'Default: <exposure folder>/catalog/. One file per '
                        'SCA (native MAST parquet), plus the combined .reg '
                        'file DS9 loads.')
    p.add_argument('--ds9-target', default=None,
                   help='DS9 XPA target name (ds9 mode). Default: first DS9 found.')
    p.add_argument('--workers', type=int, default=8,
                   help='Concurrent SCA streams (default 8). Set to 1 for '
                        'the old sequential path.')
    p.add_argument('--no-metadata', action='store_true',
                   help='Skip the per-exposure metadata CSV (written by default '
                        'alongside FITS/DS9 output, since we already streamed '
                        'the data).')
    p.add_argument('--metadata-dir', default=None,
                   help='Override the per-exposure metadata CSV location. '
                        'Default: <exposure folder>/. CSVs are dropped into '
                        'this directory directly (no per-exposure subfolder).')
    p.add_argument('--list', action='store_true',
                   help='Just list the matching exposures and exit without '
                        'streaming anything.')
    p.add_argument('--max-rows', type=int, default=50,
                   help='Max exposures to show in the summary (default 50)')

    args = p.parse_args()

    res = list_data_from_args(args)

    if res.n_exposures == 0:
        print_summary(res, max_rows=args.max_rows)
        print("\nNo exposures match — nothing to output.")
        return

    if args.list:
        print_summary(res, max_rows=args.max_rows)
        return

    indices = _resolve_exposure_keys(res, args.exposures)
    scas = None
    if args.scas is not None:
        from roman_mast import parse_int_spec
        scas = parse_int_spec(args.scas)

    write_metadata = not args.no_metadata
    _log(f"Dispatching {len(indices)} exposure(s) → {args.to}"
         + (" (+ metadata CSV)" if write_metadata else ""))

    if write_metadata:
        # Import lazily so a missing export_metadata_csv (or its deps) doesn't
        # kill FITS/DS9 output for users who don't want the CSV anyway.
        from export_metadata_csv import write_metadata_csv

    for idx in indices:
        exp = res.select(idx)

        # One folder per exposure holds every export for that exposure. Layout:
        #   <exp_dir>/
        #       metadata_{visit}_exp{NN}.csv     (unless --no-metadata)
        #       catalog/                          (ds9 mode, unless --no-catalog)
        #           r..._wfi{NN}_..._cat.parquet
        #           catalog_{visit}_exp{NN}.reg
        #       sca_{NN}.fits                     (fits mode)
        # --out-dir sets the root; each exposure gets its own v..._exp.. folder
        # under it. --metadata-dir / --catalog-dir override the sink location
        # if a user wants to route them elsewhere (they don't get the per-
        # exposure subfolder in that case — the escape hatch is meant to be
        # simple).
        root = args.out_dir or '.'
        exp_dir = os.path.join(
            root, f'v{exp.visit_id}_exp{exp.exposure:02d}',
        )
        os.makedirs(exp_dir, exist_ok=True)

        # `out_dir` is the argument to `to_fits_files` — for fits mode it's
        # exp_dir directly (SCAs land at the top of the exposure folder).
        out_dir = exp_dir if args.to == 'fits' else None

        # Where to write the metadata CSV.
        if write_metadata:
            meta_dir = args.metadata_dir if args.metadata_dir is not None else exp_dir
            os.makedirs(meta_dir, exist_ok=True)
            meta_path = os.path.join(
                meta_dir,
                f'metadata_{exp.visit_id}_exp{exp.exposure:02d}.csv',
            )

        # One stream feeds every sink. `stream_materialized` reads each SCA's
        # data arrays into memory as soon as it's opened, so sinks running
        # tens of seconds later never re-hit the 60 s pre-signed S3 URL.
        dm_dict = stream_materialized(
            exp, res.missions, scas=scas, max_workers=args.workers,
        )
        try:
            if write_metadata:
                try:
                    write_metadata_csv(dm_dict, exp, output=meta_path)
                except Exception as e:
                    _log(f"WARNING: metadata CSV failed for exp "
                         f"{exp.visit_id}/{exp.exposure:02d}: {e}")

            if args.to == 'fits':
                to_fits_files(
                    dm_dict, exp,
                    out_dir=out_dir,
                    compress=args.compress,
                    sip_degree=args.sip_degree,
                )
            else:  # ds9
                # Download L4 catalogs (parquet — can't stream, has to hit
                # disk). Skipped whole-cloth when --no-catalog is set so we
                # don't pay for the HTTP GETs users don't want. Otherwise
                # we save them into --catalog-dir (default cwd) — same
                # spirit as the metadata CSV: we're already fetching, so
                # keep the file for offline reuse.
                catalog_paths = None
                if not args.no_catalog:
                    cat_dir = (args.catalog_dir if args.catalog_dir is not None
                               else os.path.join(exp_dir, 'catalog'))
                    os.makedirs(cat_dir, exist_ok=True)
                    catalog_paths = download_catalogs(
                        exp, res.missions, scas=scas,
                        out_dir=cat_dir, max_workers=args.workers,
                    )
                    # Drop the assembled DS9 region file alongside so users
                    # can reload it without re-running the whole pipeline.
                # 'same' → single-color mode (opt out of the point/extended split).
                ext_color = (None if args.catalog_extended_color == 'same'
                             else args.catalog_extended_color)
                if any(p is not None for p in catalog_paths.values()):
                    reg_path = os.path.join(
                        cat_dir,
                        f'catalog_{exp.visit_id}_exp{exp.exposure:02d}.reg',
                    )
                    try:
                        reg_text = _build_region_text(
                            catalog_paths,
                            radius_arcsec=args.catalog_radius,
                            color=args.catalog_color,
                            extended_color=ext_color,
                            include_flagged=args.catalog_include_flagged,
                            label_mode=args.catalog_label,
                        )
                        with open(reg_path, 'w') as fh:
                            fh.write(reg_text)
                        _log(f"Wrote DS9 region file → {reg_path}")
                    except Exception as e:
                        _log(f"WARNING: could not write region file "
                             f"{reg_path}: {type(e).__name__}: {e}")
                to_ds9(
                    dm_dict, exp,
                    sip_degree=args.sip_degree,
                    dq_overlay=not args.no_dq_overlay,
                    catalog_paths=catalog_paths,
                    catalog_radius_arcsec=args.catalog_radius,
                    catalog_color=args.catalog_color,
                    catalog_extended_color=ext_color,
                    catalog_include_flagged=args.catalog_include_flagged,
                    catalog_label_mode=args.catalog_label,
                    catalog_show_labels=args.catalog_show_labels,
                    ds9_target=args.ds9_target,
                )
        finally:
            close_streams(dm_dict)


if __name__ == '__main__':
    _cli()
