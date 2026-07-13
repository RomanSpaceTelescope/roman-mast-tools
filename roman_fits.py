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

    to_ds9(af_dict, exposure, dq_overlay=True, ds9_target=None, ...)
        Pipe every SCA into a single DS9 WCS mosaic frame via XPA, optionally
        with a DQ mask overlay. No local files touched.

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


def stream_materialized(exposure: Exposure, missions, *, scas=None,
                        show_progress=True) -> dict:
    """Stream every SCA of `exposure` and pre-load its arrays inline.

    Like `roman_mast.stream_exposure`, but interleaves the S3 open with an
    immediate array read per SCA, so no consumer downstream ever needs the
    (60 s) pre-signed URL again. Use this whenever multiple sinks share the
    stream (FITS + DS9 + metadata CSV, etc.) — the plain `stream_exposure`
    variant is fine for a single-pass sink that finishes each SCA within
    the URL's lifetime.

    Returns
    -------
    dict[int, DataModel]
        Mapping SCA integer → open roman_datamodels DataModel with its
        big arrays already resident in memory. A None value marks an SCA
        we failed to stream.
    """
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

    _log(f"Streaming exposure visit_id={exposure.visit_id} "
         f"exp={exposure.exposure} ({len(pairs)} SCA files, materializing "
         f"inline)")

    dm_dict: dict = {}
    iterator = tqdm(pairs, desc="Streaming SCAs", disable=not show_progress)
    for sca, filename in iterator:
        try:
            af = missions.read_product(filename)
            dm = rdm.open(af)
            _materialize_dm(dm)
            dm_dict[sca] = dm
        except Exception as e:
            _log(f"ERROR streaming SCA {sca:02d} ({filename}): "
                 f"{type(e).__name__}: {e}")
            dm_dict[sca] = None

    return dm_dict


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
# to_ds9 — pipe an in-memory MEF via XPA
# ---------------------------------------------------------------------------

def to_ds9(
    af_dict: dict,
    exposure: Exposure,
    *,
    sip_degree: int = 4,
    dq_overlay: bool = True,
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

    d.set('zoom to fit')

    extras = f" + DQ overlay ({n_dq} SCAs)" if n_dq else ""
    _log(f"Loaded {n_data} SCAs into DS9{extras}  "
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

Examples:
  # See what's available (no output yet)
  python roman_fits.py --program 114 --pass 57 --sca-only --list

  # Write one exposure to disk as 18 FITS files
  python roman_fits.py --program 114 --pass 57 --sca-only \\
      --exposures 1 --to fits --out-dir /tmp/wfi_exp1

  # Same, RICE compressed
  python roman_fits.py --program 114 --pass 57 --sca-only \\
      --exposures 1 --to fits --compress

  # Multiple exposures — each gets its own subdirectory under --out-dir
  python roman_fits.py --program 114 --pass 57 --sca-only \\
      --exposures 1-4 --to fits --out-dir /tmp/wfi

  # Stream to DS9 (needs `ds9 &` running and pyds9 installed)
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
                   help="Base output directory (fits mode). Default: "
                        "wfi_fits_{visit}_exp{NN}/ in cwd; for multi-exposure "
                        "runs, one such folder per exposure under this base.")
    p.add_argument('--compress', action='store_true',
                   help='RICE_1 tile-compressed .fits.fz (fits mode only)')
    p.add_argument('--sip-degree', type=int, default=4,
                   help='SIP polynomial degree for gwcs → FITS (default 4)')
    p.add_argument('--no-dq-overlay', action='store_true',
                   help='Skip the DQ (bad-pixel) mask overlay in ds9 mode')
    p.add_argument('--ds9-target', default=None,
                   help='DS9 XPA target name (ds9 mode). Default: first DS9 found.')
    p.add_argument('--no-metadata', action='store_true',
                   help='Skip the per-exposure metadata CSV (written by default '
                        'alongside FITS/DS9 output, since we already streamed '
                        'the data).')
    p.add_argument('--metadata-dir', default=None,
                   help='Where to drop the metadata CSVs. Default: same as '
                        '--out-dir in fits mode, cwd in ds9 mode.')
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

    multi = len(indices) > 1
    write_metadata = not args.no_metadata
    _log(f"Dispatching {len(indices)} exposure(s) → {args.to}"
         + (" (+ metadata CSV)" if write_metadata else ""))

    if write_metadata:
        # Import lazily so a missing export_metadata_csv (or its deps) doesn't
        # kill FITS/DS9 output for users who don't want the CSV anyway.
        from export_metadata_csv import write_metadata_csv

    for idx in indices:
        exp = res.select(idx)

        # Resolve per-exposure output paths up front.
        if args.to == 'fits':
            if args.out_dir and multi:
                # Multi-exposure: one sub-folder per exposure under --out-dir.
                out_dir = os.path.join(
                    args.out_dir,
                    f'v{exp.visit_id}_exp{exp.exposure:02d}',
                )
            elif args.out_dir:
                out_dir = args.out_dir
            else:
                out_dir = f'wfi_fits_{exp.visit_id}_exp{exp.exposure:02d}'
        else:
            out_dir = None  # DS9 mode — nothing on disk unless metadata_dir set.

        # Where to write the metadata CSV.
        if write_metadata:
            if args.metadata_dir is not None:
                meta_dir = args.metadata_dir
            elif out_dir is not None:
                meta_dir = out_dir
            else:
                meta_dir = '.'
            os.makedirs(meta_dir, exist_ok=True)
            meta_path = os.path.join(
                meta_dir,
                f'metadata_{exp.visit_id}_exp{exp.exposure:02d}.csv',
            )

        # One stream feeds every sink. `stream_materialized` reads each SCA's
        # data arrays into memory as soon as it's opened, so sinks running
        # tens of seconds later never re-hit the 60 s pre-signed S3 URL.
        dm_dict = stream_materialized(exp, res.missions, scas=scas)
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
                to_ds9(
                    dm_dict, exp,
                    sip_degree=args.sip_degree,
                    dq_overlay=not args.no_dq_overlay,
                    ds9_target=args.ds9_target,
                )
        finally:
            close_streams(dm_dict)


if __name__ == '__main__':
    _cli()
