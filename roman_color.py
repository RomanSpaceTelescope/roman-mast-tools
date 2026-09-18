"""Build an RGB composite of one SCA (or the full focal plane) in DS9 and/or
as a PNG.

Each RGB channel is specified as a comma-separated key=value list of
`roman_mast.list_data` filters, plus `exposure=N` to pick which exposure
inside that query. Any filter accepted by roman_mast (program, pass,
execution_plan, segment, observation, visit, detector, optical_element,
data_level) can be used.

Example (program 1047 pass 1, three filters):

    roman-color --sca 7 \\
        --blue  program=1047,pass=1,observation=12,exposure=2 \\
        --green program=1047,pass=1,observation=5,exposure=2  \\
        --red   program=1047,pass=1,observation=9,exposure=2

Same SCA, also saved as a PNG (rgb_sca07.png in the output dir). DS9 is
optional: if it isn't running you still get the PNG.

    roman-color --sca 7 --save-png --program 1047 --pass 1 \\
        --blue observation=12,exposure=2 --green observation=5,exposure=2 \\
        --red observation=9,exposure=2

Full focal plane as a PNG with no DS9 at all. By default the SCAs are tiled
in the fixed WFI focal-plane layout (compact, no reprojection); add
--png-wcs for a north-up reprojected frame, --save-fits for stitched FITS:

    roman-color --mosaic --no-ds9 --save-png mosaic.png --program 1047 \\
        --pass 1 --blue observation=12,exposure=2 \\
        --green observation=5,exposure=2 --red observation=9,exposure=2
"""
from __future__ import annotations

import argparse
import io
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import shift as ndi_shift
from scipy.signal import fftconvolve

from roman_mast import list_data, close_streams
from roman_fits import (
    stream_materialized, mosaic_grid, stitch_mosaic, write_mosaic_fits,
    _write_hdulist,
)

# pyds9 is optional (not a declared dependency) — --no-ds9 / --save-png /
# --save-fits must work on machines without DS9. Same guard as roman_fits.
try:
    import pyds9
except ImportError:  # pragma: no cover
    pyds9 = None

DS9_TARGET = None  # None → pyds9 default; else the XPA target name


def _connect_ds9():
    if pyds9 is None:
        raise ImportError(
            "pyds9 is required to push into DS9. Install with `pip install "
            "pyds9` (and run `ds9 &`), or pass --no-ds9 with --save-png / "
            "--save-fits for file output only."
        )
    print('[rgb] connecting to DS9 via pyds9', file=sys.stderr)
    return pyds9.DS9(target=DS9_TARGET) if DS9_TARGET else pyds9.DS9()

# Preferred output layout: mounted shared storage on server nodes. Falls
# back to CWD-relative if the mount isn't present (e.g. laptop dev).
_PREFERRED_OUT_MOUNT = '/mnt/roman-science-east-2'
_PREFERRED_OUT_SUBDIR = 'mrizzo'


def _default_out_dir(subdir: str) -> str:
    """Return an absolute path under the preferred output root if the shared
    mount is present, else under the current working directory. Creates the
    per-user subdir and any intermediate program/leaf dirs as needed."""
    if os.path.isdir(_PREFERRED_OUT_MOUNT):
        base = os.path.join(_PREFERRED_OUT_MOUNT, _PREFERRED_OUT_SUBDIR,
                            subdir)
    else:
        base = os.path.abspath(subdir)
    os.makedirs(base, exist_ok=True)
    return base


def _spec_dir_name(specs: dict, *, sca: int | None = None) -> str:
    """Build a descriptive folder name from the three channel specs.

    Format: p{program}_pass{pass}_exp{exp}_B{filter}-G{filter}-R{filter}[_sca{NN}]

    Fields the specs don't share are collapsed to a placeholder rather
    than duplicating. Filter names are pulled from `filter=` /
    `optical_element=` in each spec.
    """
    def _flt(spec):
        return str(spec.get('optical_element',
                            spec.get('filter', 'UNK'))).upper()

    def _shared(key):
        vals = {spec.get(key) for spec in specs.values() if key in spec}
        if len(vals) == 1:
            v = next(iter(vals))
            return v if v is not None else None
        return None  # differ across channels — omit

    parts = []
    pass_ = _shared('pass_')
    if pass_ is not None:
        parts.append(f'pass{int(pass_):03d}')
    expn = _shared('exposure')
    if expn is not None:
        parts.append(f'exp{int(expn):02d}')

    filters = (f"B{_flt(specs['blue'])}-"
               f"G{_flt(specs['green'])}-"
               f"R{_flt(specs['red'])}")
    parts.append(filters)

    if sca is not None:
        parts.append(f'sca{sca:02d}')

    leaf = '_'.join(parts) if parts else 'rgb'

    # Nest under a program-level parent folder when program is shared.
    prog = _shared('program')
    if prog is not None:
        return os.path.join(f'{int(prog):05d}', leaf)
    return leaf

# roman_mast.list_data filter keys we accept in --red/--green/--blue specs.
# `pass` is a Python keyword → mapped to `pass_` when passed to list_data.
FILTER_KEYS = {
    'program', 'pass', 'pass_', 'execution_plan', 'segment', 'observation',
    'visit', 'detector', 'optical_element', 'data_level',
}
INT_KEYS = {'program', 'pass', 'pass_', 'execution_plan', 'segment',
            'observation', 'visit', 'exposure'}


def parse_spec(spec: str) -> dict:
    """Parse `key=val,key=val,...` into a dict, coercing ints where sensible.

    Accepts `filter=Fxxx` as a shortcut for `optical_element=Fxxx`. When
    filter is set and observation is not, roman-color picks the first
    exposure returned by list_data that matches the filter + exposure.
    """
    out = {}
    for part in spec.split(','):
        part = part.strip()
        if not part:
            continue
        if '=' not in part:
            raise ValueError(f'Bad spec fragment {part!r}: expected key=value')
        k, v = part.split('=', 1)
        k = k.strip()
        v = v.strip()
        if k == 'pass':
            k = 'pass_'
        if k == 'filter':
            k = 'optical_element'
        if k in INT_KEYS or (k == 'pass_'):
            v = int(v)
        out[k] = v
    if 'exposure' not in out:
        raise ValueError(f'Spec {spec!r} needs exposure=N')
    return out


def _resolve_exposure(spec: dict):
    """Query MAST and return (exposure, missions, filters_used)."""
    exposure_number = int(spec['exposure'])
    filters = {k: v for k, v in spec.items()
               if k != 'exposure' and (k in FILTER_KEYS or k == 'pass_')}

    res = list_data(**filters)
    if res.n_exposures == 0:
        raise RuntimeError(f'No exposures returned for {filters}')

    matches = [c for c in res.exposures if int(c.exposure) == exposure_number]
    if not matches:
        avail = sorted({int(e.exposure) for e in res.exposures})
        raise RuntimeError(
            f'{filters}: no exposure {exposure_number}; available: {avail}'
        )
    if len(matches) > 1:
        print(f'[rgb] WARNING: {filters} exposure={exposure_number} matches '
              f'{len(matches)} exposures; picking the first '
              f'(observation={getattr(matches[0], "observation", "?")})',
              file=sys.stderr)
    exp = matches[0]
    print(f'[rgb] {filters}  exposure={exposure_number}  '
          f'filter={exp.optical_element}  visit_id={exp.visit_id}',
          file=sys.stderr)
    return exp, res.missions, filters


def _log_filenames(exp, scas, tag: str):
    """Print the exact _cal.asdf filenames roman-color is about to stream."""
    wanted = set(scas) if scas is not None else set(exp.scas)
    print(f'[rgb] {tag}: {len(wanted)} file(s) to stream:', file=sys.stderr)
    for sca, fname in zip(exp.scas, exp.filenames):
        if sca in wanted:
            print(f'[rgb]   SCA{sca:02d}: {fname}', file=sys.stderr)


def fetch_sca(spec: dict, sca: int):
    """Return (data, wcs_header, filter_name, obs, expn, visit_id)."""
    exp, missions, _ = _resolve_exposure(spec)
    exposure_number = int(spec['exposure'])
    _log_filenames(exp, [sca],
                   tag=f'exp={exposure_number} filter={exp.optical_element}')
    dm_dict = stream_materialized(exp, missions, scas=[sca], max_workers=1)
    dm = dm_dict.get(sca)
    if dm is None:
        raise RuntimeError(f'SCA {sca} missing for exp={exposure_number}')
    try:
        data = np.asarray(dm.data[...], dtype=np.float32)
        gwcs = dm.meta.wcs
        wcs_hdr = gwcs.to_fits_sip(bounding_box=gwcs.bounding_box, degree=4)
        filt = str(exp.optical_element)
    finally:
        close_streams(dm_dict)
    obs = int(getattr(exp, 'observation', 0) or 0)
    return data, wcs_hdr, filt, obs, exposure_number, exp.visit_id


def fetch_exposure(spec: dict, scas, max_workers: int = 8):
    """Stream all requested SCAs of one exposure in parallel.

    Returns a dict keyed by SCA number:
        {sca: {'data': ndarray, 'hdr': fits.Header}}
    plus a metadata dict with filter/obs/exp/visit_id (identical across SCAs).
    """
    exp, missions, _ = _resolve_exposure(spec)
    _log_filenames(exp, scas,
                   tag=f'exp={int(exp.exposure)} '
                       f'filter={exp.optical_element}')
    dm_dict = stream_materialized(exp, missions, scas=scas,
                                  max_workers=max_workers)
    per_sca = {}
    try:
        for sca, dm in dm_dict.items():
            if dm is None:
                print(f'[rgb] WARNING: SCA {sca} stream returned None',
                      file=sys.stderr)
                continue
            data = np.asarray(dm.data[...], dtype=np.float32)
            gwcs = dm.meta.wcs
            hdr = gwcs.to_fits_sip(bounding_box=gwcs.bounding_box, degree=4)
            per_sca[sca] = {'data': data, 'hdr': hdr}
    finally:
        close_streams(dm_dict)
    meta = {
        'filter': str(exp.optical_element),
        'obs': int(getattr(exp, 'observation', 0) or 0),
        'exp': int(exp.exposure),
        'visit_id': exp.visit_id,
    }
    return per_sca, meta


def cache_name(channel: str, meta: dict, sca: int) -> str:
    return (f'{channel}_obs{meta["obs"]:02d}_exp{meta["exp"]:02d}_'
            f'{meta["filter"]}_sca{sca:02d}.fits')


def write_fits(path: str, data: np.ndarray, wcs_hdr: fits.Header):
    hdu = fits.PrimaryHDU(data=data.astype(np.float32), header=wcs_hdr)
    _write_hdulist(fits.HDUList([hdu]), path, overwrite=True)
    print(f'[rgb] wrote {path}', file=sys.stderr)


def _hdulist_bytes(hdulist: fits.HDUList) -> bytes:
    """Serialize an HDUList to bytes.

    writeto(BytesIO) never touches a file descriptor, so it cannot raise the
    'Bad file descriptor' error that on-disk writeto() hits on Mountpoint-S3
    mounts (see roman_fits._write_hdulist). The tempfile fallback is kept as
    a belt-and-braces path only.
    """
    import tempfile
    try:
        buf = io.BytesIO()
        hdulist.writeto(buf)
        return buf.getvalue()
    except OSError:
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as tmp:
            tmp_path = tmp.name
        try:
            hdulist.writeto(tmp_path, overwrite=True)
            with open(tmp_path, 'rb') as f:
                return f.read()
        finally:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass


def _fits_bytes(data: np.ndarray, hdr: fits.Header) -> bytes:
    hdulist = fits.HDUList([fits.PrimaryHDU(data=data.astype(np.float32),
                                            header=hdr)])
    return _hdulist_bytes(hdulist)


def _clean_for_xcorr(img, saturate_sigma: float = 10.0):
    """Prep an image for FFT cross-correlation.

    Subtracts the sigma-clipped median, divides by the sigma-clipped std
    (so each channel contributes similar weight regardless of intensity
    scale), then clips at ±`saturate_sigma`. The clip is critical: on
    fields with a few very bright saturated stars, the raw correlation
    peak parks at (0, 0) because saturated pixels dominate the sum via
    self-match, drowning out the geometric-shift signal from the fainter
    star field. Clipping to a few sigma flattens the bright peaks so the
    correlation is driven by the many fainter sources whose positions
    actually shift between exposures.
    """
    arr = np.array(img, dtype=np.float32)
    finite = np.isfinite(arr)
    if finite.any():
        _, med, std = sigma_clipped_stats(arr[finite])
    else:
        med, std = 0.0, 1.0
    std = max(std, 1e-6)
    out = np.where(finite, (arr - med) / std, 0.0)
    if saturate_sigma:
        out = np.clip(out, -saturate_sigma, saturate_sigma)
    return out.astype(np.float32)


SEARCH_RADIUS_PX = 10


def _bounded_shift(ref_full, mov_full, search_radius=SEARCH_RADIUS_PX):
    """Bounded FFT cross-correlation. Returns (dy, dx, peak)."""
    ref = _clean_for_xcorr(ref_full)
    mov = _clean_for_xcorr(mov_full)
    corr = fftconvolve(ref, mov[::-1, ::-1], mode='same')
    cy, cx = corr.shape[0] // 2, corr.shape[1] // 2
    y0, y1 = cy - search_radius, cy + search_radius + 1
    x0, x1 = cx - search_radius, cx + search_radius + 1
    window = corr[y0:y1, x0:x1]
    py, px = np.unravel_index(np.argmax(window), window.shape)
    dy = (py + y0) - cy
    dx = (px + x0) - cx
    return float(dy), float(dx), float(window.max())


def align_to_ref(raw: np.ndarray, ref_aligned: np.ndarray, *,
                 label: str = ''):
    """Shift `raw` to match `ref_aligned` using bounded FFT correlation."""
    finite_in = np.isfinite(raw)
    dy, dx, peak = _bounded_shift(ref_aligned, raw)
    if label:
        print(f'[rgb] {label} shift (dy, dx) = ({dy:+.1f}, {dx:+.1f}) px  '
              f'peak={peak:.3g}', file=sys.stderr)
    clean = np.where(finite_in, raw, 0.0)
    shifted = ndi_shift(clean, shift=(dy, dx), order=3,
                        mode='constant', cval=np.nan)
    finite_shifted = ndi_shift(
        finite_in.astype(np.float32), shift=(dy, dx), order=0,
        mode='constant', cval=0.0,
    ) > 0.5
    shifted[~finite_shifted] = np.nan
    return shifted.astype(np.float32)


# Per-channel asinh stretch: limits = [median - LOW_SIGMA*rms,
# median + GAIN_SIGMA[channel]*rms]. Red (F184) is pushed harder to make up
# for its lower throughput so the composite doesn't go blue-dominated.
# Shared by the DS9 pushes and the PNG renders in both modes.
GAIN_SIGMA = {'blue': 30.0, 'green': 25.0, 'red': 20.0}
LOW_SIGMA = 1.0


def _sca_limits(arr, channel: str, *, tag: str = '') -> tuple[float, float]:
    """(lo, hi) asinh limits for one channel of one SCA from sigma-clipped
    stats of its finite pixels. Logs the numbers so DS9 and PNG runs are
    comparable."""
    finite = np.isfinite(arr)
    if finite.any():
        _, med, std = sigma_clipped_stats(arr[finite])
    else:
        med, std = 0.0, 1.0
    lo = float(med - LOW_SIGMA * std)
    hi = float(med + GAIN_SIGMA[channel] * std)
    print(f'[rgb] {tag}{channel}: bkg={med:.4g}  rms={std:.4g}  '
          f'limits=[{lo:.4g}, {hi:.4g}]  asinh', file=sys.stderr)
    return lo, hi


def _run_ds9(layers, blue):
    """Push the three aligned layers into DS9 as an RGB frame."""
    d = _connect_ds9()
    d.set('frame delete all')
    d.set('frame new rgb')
    for channel in ('red', 'green', 'blue'):
        lyr = layers[channel]
        d.set(f'rgb channel {channel}')
        d.set('fits', _fits_bytes(lyr['aligned'], blue['hdr']))
        lo, hi = _sca_limits(lyr['aligned'], channel, tag='ds9 ')
        d.set(f'scale limits {lo} {hi}')
        d.set('scale asinh')
    d.set('rgb channel red')
    d.set('wcs align no')  # image-aligned display; WCS not yet trustworthy
    d.set('zoom to fit')


def _sca_png(layers, png_path: str):
    """Render the three aligned single-SCA layers to an RGB PNG.

    Same per-channel asinh limits as the DS9 push (`_sca_limits`), array
    row 0 at the bottom like DS9's default display and the mosaic PNGs.
    Needs only `layers[channel]['aligned']`, so it works for both the
    streamed and the --from-cache layer dicts.
    """
    import matplotlib.pyplot as plt
    rgb = None
    for idx, channel in enumerate(('red', 'green', 'blue')):
        arr = layers[channel]['aligned']
        lo, hi = _sca_limits(arr, channel, tag='png ')
        if rgb is None:
            rgb = np.zeros(arr.shape + (3,), dtype=np.float32)
        rgb[:, :, idx] = _asinh_scale(arr, lo, hi)
    rgb = np.flipud(rgb)
    fig, ax = plt.subplots(figsize=(12, 12), dpi=150)
    ax.imshow(rgb, origin='upper', interpolation='nearest')
    ax.set_axis_off()
    fig.savefig(png_path, bbox_inches='tight', pad_inches=0,
                facecolor='black')
    plt.close(fig)
    print(f'[rgb] wrote PNG {png_path}', file=sys.stderr)


def _push_ds9(push, *args, file_outputs: bool, **kwargs) -> bool:
    """Run a DS9 push (`_run_ds9` / `_run_ds9_mosaic`).

    A missing pyds9, no running DS9, or no display all surface here as an
    exception. When the user asked for file products (--save-png /
    --save-fits) those are already on disk by the time we get here, so
    the run is a success from the file point of view: log a WARNING and
    return False. With no file outputs requested the exception propagates,
    because the user asked for DS9 and got nothing.
    """
    try:
        push(*args, **kwargs)
    except Exception as e:  # pyds9 raises assorted types on XPA failure
        if not file_outputs:
            raise
        print(f'[rgb] WARNING: DS9 push failed ({type(e).__name__}: {e}); '
              'file outputs were written, continuing.', file=sys.stderr)
        return False
    return True


ALL_SCAS = list(range(1, 19))


def _stretch_per_channel(all_arrays_by_channel):
    """Compute one (lo, hi) per channel from pooled sigma-clipped stats.

    all_arrays_by_channel : {channel: [aligned_ndarray, ...]}
        Concatenated across SCAs, so the mosaic gets a uniform stretch.
        Uses the module-level GAIN_SIGMA / LOW_SIGMA.
    """
    limits = {}
    for channel, arrays in all_arrays_by_channel.items():
        # Sample from every SCA — full stack would be huge. Use every 8th px.
        samples = []
        for arr in arrays:
            f = np.isfinite(arr)
            if f.any():
                samples.append(arr[f][::8])
        if samples:
            pool = np.concatenate(samples)
            _, med, std = sigma_clipped_stats(pool)
        else:
            med, std = 0.0, 1.0
        limits[channel] = (float(med - LOW_SIGMA * std),
                           float(med + GAIN_SIGMA[channel] * std))
        print(f'[rgb] mosaic {channel}: bkg={med:.4g}  rms={std:.4g}  '
              f'limits={limits[channel]}', file=sys.stderr)
    return limits


def _build_channel_mef(sca_layers, ref_hdrs, channel: str,
                       out_path: str | None = None) -> bytes:
    """Build a multi-extension FITS (one ImageHDU per SCA) for one RGB channel.

    Every SCA extension carries its own SIP WCS in the header, so DS9's
    `fits mosaicimage wcs` command stitches them into a single WCS-aligned
    view. Returns the serialized MEF as bytes; also writes to disk if
    `out_path` is provided.
    """
    hdul = fits.HDUList([fits.PrimaryHDU()])
    for sca in sorted(sca_layers):
        arr = sca_layers[sca].get(channel)
        if arr is None:
            continue
        hdr = ref_hdrs[sca].copy()
        hdr['EXTNAME'] = f'SCA{sca:02d}'
        hdr['SCANUM'] = (sca, 'SCA number (1-18)')
        hdr['CHANNEL'] = (channel, 'RGB channel')
        hdul.append(fits.ImageHDU(data=arr.astype(np.float32),
                                  header=hdr, name=f'SCA{sca:02d}'))
    # Serialize once in memory; when we have a target path push the bytes out
    # as a single sequential write (the "Bad file descriptor" seen on the RES
    # S3 mounts comes from astropy's on-disk writeto dup()/close()-ing the fd
    # mid-file, which completes the Mountpoint-S3 upload early — see
    # roman_fits._write_hdulist). The same bytes feed the DS9 pipe.
    if out_path is not None:
        _write_hdulist(hdul, out_path, overwrite=True)
        with open(out_path, 'rb') as f:
            data = f.read()
        print(f'[rgb] wrote {out_path} ({len(data)/1e6:.0f} MB, '
              f'{len(hdul)-1} SCAs)', file=sys.stderr)
    else:
        data = _hdulist_bytes(hdul)
    return data


def _asinh_scale(arr, lo, hi):
    """Match DS9's asinh: normalize to [0, 1] then apply the asinh stretch."""
    x = (arr - lo) / max(hi - lo, 1e-12)
    x = np.clip(np.nan_to_num(x, nan=0.0), 0.0, 1.0)
    return np.arcsinh(10.0 * x) / np.arcsinh(10.0)


def _focal_plane_png(sca_layers, limits, png_path, *, decimate=4):
    """Render the RGB mosaic in the fixed WFI focal-plane layout (no WCS).

    Each SCA is block-averaged by `decimate`, stretched with the shared
    per-channel asinh limits, and dropped onto a black canvas at its
    focal-plane position from `roman_phot._WFI_SCA_LAYOUT`. Nothing is
    reprojected: the canvas is the physical detector footprint (~266 x
    165 mm), so the frame has no blank corners from the roll angle, is the
    same orientation every time (+X focal-plane to the left, +Y up, same
    convention as roman_phot's mosaics), and peak memory is one decimated
    tile rather than three full-focal-plane arrays.
    """
    from roman_phot import (_WFI_SCA_LAYOUT, _ROMAN_SCA_FULL_SIZE,
                            _ROMAN_PIXEL_SCALE_MM)
    import matplotlib.pyplot as plt

    mm_per_px = _ROMAN_PIXEL_SCALE_MM * decimate        # canvas pixel pitch
    half = _ROMAN_SCA_FULL_SIZE * _ROMAN_PIXEL_SCALE_MM / 2  # 20.48 mm
    xs = [cx for cx, _, _ in _WFI_SCA_LAYOUT.values()]
    ys = [cy for _, cy, _ in _WFI_SCA_LAYOUT.values()]
    x_hi, y_hi = max(xs) + half, max(ys) + half
    x_lo, y_lo = min(xs) - half, min(ys) - half
    width = int(np.ceil((x_hi - x_lo) / mm_per_px))
    height = int(np.ceil((y_hi - y_lo) / mm_per_px))
    canvas = np.zeros((height, width, 3), dtype=np.uint8)

    def _block_mean(arr, d):
        h, w = (arr.shape[0] // d) * d, (arr.shape[1] // d) * d
        blocks = arr[:h, :w].reshape(h // d, d, w // d, d)
        with warnings.catch_warnings():
            # All-NaN blocks (DQ holes) become NaN -> black; no warning.
            warnings.simplefilter('ignore', RuntimeWarning)
            return np.nanmean(blocks, axis=(1, 3))

    n_tiles = 0
    for sca in sorted(sca_layers):
        cx, cy, _ = _WFI_SCA_LAYOUT[sca]
        # Canvas column grows toward -X (East-left), row grows toward -Y.
        # Raw array row 0 lands at the bottom of the tile: the same flipud
        # the WCS path applies to the whole mosaic.
        col0 = int(round((x_hi - (cx + half)) / mm_per_px))
        row0 = int(round((y_hi - (cy + half)) / mm_per_px))
        for idx, channel in enumerate(('red', 'green', 'blue')):
            arr = sca_layers[sca].get(channel)
            if arr is None:
                continue
            lo, hi = limits[channel]
            tile = np.flipud(_asinh_scale(_block_mean(arr, decimate), lo, hi))
            h, w = tile.shape
            r1, c1 = min(row0 + h, height), min(col0 + w, width)
            canvas[row0:r1, col0:c1, idx] = (
                255 * tile[:r1 - row0, :c1 - col0]).astype(np.uint8)
        n_tiles += 1

    fig_w = 16.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_w * height / width), dpi=150)
    ax.imshow(canvas, origin='upper', interpolation='nearest')
    ax.set_axis_off()
    fig.savefig(png_path, bbox_inches='tight', pad_inches=0,
                facecolor='black')
    plt.close(fig)
    print(f'[rgb] wrote PNG {png_path} (focal-plane layout, '
          f'{n_tiles} SCAs, canvas {width}x{height} px at decimate='
          f'{decimate})', file=sys.stderr)


def _stitch_channels(sca_layers, limits, ref_hdrs, *, fits_dir=None,
                     png_path=None, channel_meta=None, png_wcs=False):
    """Stitch each channel's per-SCA aligned arrays into one focal-plane image.

    All three channels land on the same celestial grid (derived from blue's
    per-SCA WCS headers). Per channel, optionally:
      - write the stitched image as plain float32 FITS to
        `fits_dir/stitched_{channel}.fits` (uncompressed — the DS9 MEFs
        written by _run_ds9_mosaic are per-SCA tiles, not a single image);
      - fold it into an RGB PNG (same asinh + per-channel limits as the DS9
        push) saved to `png_path` via matplotlib, no DS9 needed. Only when
        `png_wcs` is set; the default PNG is the fixed focal-plane layout
        from `_focal_plane_png`, which needs no reprojection.

    Channels are processed one at a time and dropped before the next so
    only one ~0.8 GB full-focal-plane array is resident.
    """
    if png_path is not None and not png_wcs:
        # Default PNG: fixed focal-plane tiling. Cheap, and no blank corners
        # from putting a rolled focal plane on a north-up sky grid.
        _focal_plane_png(sca_layers, limits, png_path)
        png_path = None
    if fits_dir is None and png_path is None:
        return

    def _tiles(channel):
        return [(sca_layers[sca][channel], ref_hdrs[sca])
                for sca in sorted(sca_layers)
                if channel in sca_layers[sca]]

    # Common grid from blue's tiles; green/red are pixel-aligned to blue
    # so they share its WCS.
    common_wcs = common_shape = None
    rgb = None
    for idx, channel in enumerate(('red', 'green', 'blue')):
        tiles = _tiles(channel)
        if not tiles:
            print(f'[rgb] stitch: no {channel} tiles; skipping',
                  file=sys.stderr)
            continue
        if common_wcs is None:
            common_wcs, common_shape = mosaic_grid(_tiles('blue') or tiles)
        print(f'[rgb] stitching {channel} ({len(tiles)} SCAs)',
              file=sys.stderr)
        mosaic, _, _, _ = stitch_mosaic(
            tiles, common_wcs=common_wcs, shape_out=common_shape,
        )

        if fits_dir is not None:
            extra = {
                'CHANNEL': (channel, 'RGB channel'),
                'NTILES': (len(tiles), 'Number of SCA tiles stitched'),
                'COMBINE': ('mean', 'reproject_and_coadd combine_function'),
            }
            meta = (channel_meta or {}).get(channel)
            if meta:
                extra['FILTER'] = (meta['filter'], 'Optical element / filter')
                extra['VISITID'] = (meta['visit_id'], 'Roman visit ID')
                extra['EXPNUM'] = (meta['exp'], 'Exposure number within visit')
                extra['OBSNUM'] = (meta['obs'], 'Observation number')
            write_mosaic_fits(
                os.path.join(fits_dir, f'stitched_{channel}.fits'),
                mosaic, common_wcs, extra=extra,
            )

        if png_path is not None:
            if rgb is None:
                rgb = np.zeros((common_shape[0], common_shape[1], 3),
                               dtype=np.float32)
            lo, hi = limits[channel]
            rgb[:, :, idx] = _asinh_scale(mosaic, lo, hi)
        del mosaic

    if png_path is not None and rgb is not None:
        import matplotlib.pyplot as plt
        # Matplotlib expects sky-north-up. Astropy WCS puts DEC increasing
        # in array Y; flip vertically so the PNG looks like DS9's default.
        rgb = np.flipud(rgb)
        fig, ax = plt.subplots(figsize=(12, 12), dpi=150)
        ax.imshow(rgb, origin='upper')
        ax.set_axis_off()
        fig.savefig(png_path, bbox_inches='tight', pad_inches=0,
                    facecolor='black')
        plt.close(fig)
        print(f'[rgb] wrote PNG {png_path} (WCS grid)', file=sys.stderr)


def _run_ds9_mosaic(sca_layers, limits, ref_hdrs, out_dir=None):
    """Push a single WCS-stitched RGB mosaic into DS9.

    Builds three MEFs (one per channel, each holding all 18 SCAs as
    extensions with SIP WCS headers) and pipes each into DS9's
    `fits mosaicimage wcs` command on the matching rgb channel of one
    RGB frame. Every channel of the frame becomes a full-focal-plane
    mosaic, and DS9 combines the three into a single color view.

    Optionally also writes the three MEFs to `out_dir/mosaic_{rgb}.fits`
    so downstream tools can inspect them.
    """
    d = _connect_ds9()
    d.set('frame delete all')
    d.set('frame new rgb')

    for channel in ('red', 'green', 'blue'):
        out_path = (os.path.join(out_dir, f'mosaic_{channel}.fits')
                    if out_dir else None)
        mef = _build_channel_mef(sca_layers, ref_hdrs, channel,
                                 out_path=out_path)
        print(f'[rgb] piping {channel} mosaic ({len(mef)/1e6:.0f} MB) into DS9',
              file=sys.stderr)
        d.set(f'rgb channel {channel}')
        d.set('fits mosaicimage wcs', mef)
        lo, hi = limits[channel]
        d.set(f'scale limits {lo} {hi}')
        d.set('scale asinh')

    d.set('rgb channel red')
    # Default alignment: image coords, not WCS. Roman commissioning WCS
    # isn't calibrated yet, so image-aligned display matches the raw pixel
    # frames we already aligned by FFT.
    d.set('wcs align no')
    d.set('zoom to fit')


def run_mosaic(specs, out_dir, *, workers=8, from_cache=False,
               save_png=None, save_fits=False, no_ds9=False, png_wcs=False):
    """Stream all 18 SCAs for each of the three channels, align per-SCA
    (blue as pivot), then: write per-SCA aligned FITS, optionally stitch
    each channel into one full-focal-plane image (`save_fits` → plain FITS,
    `save_png` → RGB PNG in the fixed focal-plane layout, or on the common
    WCS grid when `png_wcs`), and push one WCS-stitched RGB frame into DS9
    unless `no_ds9`. File products are written before the DS9 push, and a
    DS9 failure after them is a warning, not an error."""
    import glob as _glob

    if from_cache:
        # Load 18×3 aligned FITS from disk.
        sca_layers = {}
        ref_hdrs = {}
        pool = {'red': [], 'green': [], 'blue': []}
        for sca in ALL_SCAS:
            layers_sca = {}
            hdr_ref = None
            for channel in ('red', 'green', 'blue'):
                pattern = os.path.join(
                    out_dir, f'{channel}_*_sca{sca:02d}.fits',
                )
                matches = sorted(_glob.glob(pattern))
                if not matches:
                    print(f'[rgb] mosaic --from-cache: no {channel} FITS for '
                          f'SCA {sca:02d}', file=sys.stderr)
                    continue
                with fits.open(matches[0]) as hdul:
                    layers_sca[channel] = np.asarray(hdul[0].data,
                                                     dtype=np.float32)
                    hdr = hdul[0].header
                if channel == 'blue':
                    hdr_ref = hdr
                pool[channel].append(layers_sca[channel])
            if layers_sca:
                sca_layers[sca] = layers_sca
                if hdr_ref is not None:
                    ref_hdrs[sca] = hdr_ref
        if not sca_layers:
            raise FileNotFoundError(
                f'--from-cache --mosaic: no aligned FITS found in {out_dir}'
            )
        limits = _stretch_per_channel(pool)
        _stitch_channels(sca_layers, limits, ref_hdrs,
                         fits_dir=out_dir if save_fits else None,
                         png_path=save_png, png_wcs=png_wcs)
        if not no_ds9:
            _push_ds9(_run_ds9_mosaic, sca_layers, limits, ref_hdrs,
                      out_dir=out_dir,
                      file_outputs=bool(save_png or save_fits))
        print(f'[rgb] done (mosaic cache).', file=sys.stderr)
        return

    # Stream: one call per filter, each pulling all 18 SCAs in parallel.
    channel_data = {}
    channel_meta = {}
    for channel, spec in specs.items():
        print(f'[rgb] mosaic: streaming {channel} '
              f'(18 SCAs, {workers} workers)', file=sys.stderr)
        per_sca, meta = fetch_exposure(spec, ALL_SCAS, max_workers=workers)
        channel_data[channel] = per_sca
        channel_meta[channel] = meta

    # Per-SCA alignment (blue is pivot). Parallelized: each SCA's
    # (green→blue, red→blue) alignment is independent, and scipy's FFT
    # releases the GIL during the transform, so threads work here without
    # the pickling overhead of a ProcessPoolExecutor.
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def _align_one_sca(sca):
        blue = channel_data['blue'].get(sca)
        if blue is None:
            return sca, None, None
        blue_arr = blue['data'].astype(np.float32)
        out = {'blue': blue_arr}
        for channel in ('green', 'red'):
            src = channel_data[channel].get(sca)
            if src is None:
                continue
            out[channel] = align_to_ref(src['data'], blue_arr,
                                        label=f'SCA{sca:02d} {channel}')
        return sca, out, blue['hdr']

    sca_layers = {}
    ref_hdrs = {}
    pool = {'red': [], 'green': [], 'blue': []}

    align_workers = min(workers, len(ALL_SCAS))
    print(f'[rgb] mosaic: aligning 18 SCAs in parallel '
          f'(workers={align_workers})', file=sys.stderr)
    with ThreadPoolExecutor(max_workers=align_workers) as ex:
        futures = [ex.submit(_align_one_sca, sca) for sca in ALL_SCAS]
        for fut in as_completed(futures):
            sca, layers_sca, hdr = fut.result()
            if layers_sca is None:
                print(f'[rgb] mosaic: SCA {sca:02d} missing in blue; skipping',
                      file=sys.stderr)
                continue
            sca_layers[sca] = layers_sca
            ref_hdrs[sca] = hdr
            for channel, arr in layers_sca.items():
                pool[channel].append(arr)

    # Write aligned FITS to disk (serial — I/O bound and cheap next to align).
    for sca in sorted(sca_layers):
        layers_sca = sca_layers[sca]
        hdr = ref_hdrs[sca]
        for channel, arr in layers_sca.items():
            name = cache_name(channel, {
                'obs': channel_meta[channel]['obs'],
                'exp': channel_meta[channel]['exp'],
                'filter': channel_meta[channel]['filter'],
            }, sca)
            write_fits(os.path.join(out_dir, name), arr, hdr)

    # One stretch shared across all SCAs so the mosaic looks uniform.
    limits = _stretch_per_channel(pool)
    _stitch_channels(sca_layers, limits, ref_hdrs,
                     fits_dir=out_dir if save_fits else None,
                     png_path=save_png, channel_meta=channel_meta,
                     png_wcs=png_wcs)
    if not no_ds9:
        _push_ds9(_run_ds9_mosaic, sca_layers, limits, ref_hdrs,
                  out_dir=out_dir, file_outputs=bool(save_png or save_fits))
    print(f'[rgb] mosaic done. FITS in {out_dir}', file=sys.stderr)


# Sentinel for `--save-png` given without a path: main() substitutes a
# mode-specific default filename inside --out-dir.
_AUTO_PNG = '__auto__'


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--sca', type=int, default=None,
                    help='SCA number (1-18). Required unless --mosaic is set.')
    ap.add_argument('--mosaic', action='store_true',
                    help='Build the full 18-SCA color mosaic. Streams every '
                         'SCA in each of the three filters, aligns each SCA '
                         'independently (blue as pivot), writes per-SCA '
                         'aligned FITS, and loads one WCS-stitched RGB frame '
                         'into DS9 (unless --no-ds9). Add --save-png and/or '
                         '--save-fits for stitched file output.')
    ap.add_argument('--workers', type=int, default=8,
                    help='Concurrent SCA streams per filter (default 8).')
    ap.add_argument('--program', type=int, default=None,
                    help='Global program ID inherited by all three layers '
                         '(override per-layer via the spec).')
    ap.add_argument('--pass', dest='pass_', type=int, default=None,
                    help='Global pass number inherited by all three layers.')
    ap.add_argument('--execution-plan', type=int, default=None,
                    help='Global execution-plan number.')
    ap.add_argument('--segment', type=int, default=None,
                    help='Global segment number.')
    ap.add_argument('--observation', type=int, default=None,
                    help='Global observation number inherited by all layers.')
    ap.add_argument('--visit', type=int, default=None,
                    help='Global visit number inherited by all layers.')
    ap.add_argument('--red', required=True,
                    help='Red-layer spec, e.g. '
                         'observation=9,exposure=2  or  filter=F184,exposure=2')
    ap.add_argument('--green', required=True, help='Green-layer spec')
    ap.add_argument('--blue', required=True, help='Blue-layer spec')
    ap.add_argument('--out-dir', default=None,
                    help='Output directory for aligned FITS and PNGs. '
                         'Default: a descriptive p<program>/<pass_exp_'
                         'filters[_scaNN|_mosaic]> folder under the shared '
                         'science mount when it exists, else under the '
                         'current directory.')
    ap.add_argument('--from-cache', action='store_true',
                    help='Skip MAST streaming + alignment; load the aligned '
                         'FITS files from --out-dir, then save/push as '
                         'requested.')
    ap.add_argument('--save-png', nargs='?', const=_AUTO_PNG, default=None,
                    metavar='PATH',
                    help='Write an RGB PNG with matplotlib using the same '
                         'asinh stretch as the DS9 view. Works with or '
                         'without DS9, in both --sca and --mosaic modes. '
                         'PATH is absolute or relative to --out-dir; with '
                         'no PATH it is rgb_sca<NN>.png / rgb_mosaic.png in '
                         '--out-dir. With --mosaic the 18 SCAs are tiled in '
                         'the fixed WFI focal-plane layout (no reprojection) '
                         'unless --png-wcs is given.')
    ap.add_argument('--png-wcs', action='store_true',
                    help='(--mosaic only) With --save-png, reproject the '
                         'PNG onto a common north-up celestial grid instead '
                         'of the fixed focal-plane layout. Holds three '
                         'full-focal-plane arrays in RAM and pads the frame '
                         'with blank corners from the roll angle.')
    ap.add_argument('--save-fits', action='store_true',
                    help='(--mosaic only) Stitch each channel\'s 18 aligned '
                         'SCAs into one full-focal-plane image on a common '
                         'WCS and write it as plain uncompressed float32 '
                         'FITS: stitched_{red,green,blue}.fits in --out-dir.')
    ap.add_argument('--no-ds9', action='store_true',
                    help='Skip the DS9 push entirely (both modes). Not '
                         'required just to get files: if DS9 is unavailable '
                         'after --save-png / --save-fits were written, the '
                         'run warns and exits cleanly anyway.')
    args = ap.parse_args()

    def _resolve_png_path(p, base, auto_name):
        """None -> None; the --save-png sentinel -> auto_name in base;
        relative -> under base; absolute unchanged. Creates the parent."""
        if p is None:
            return None
        if p == _AUTO_PNG:
            p = auto_name
        path = p if os.path.isabs(p) else os.path.join(base, p)
        os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        return path

    if not args.mosaic and args.sca is None:
        ap.error('either --sca N or --mosaic is required')

    globals_ = {}
    for k, v in (('program', args.program), ('pass_', args.pass_),
                 ('execution_plan', args.execution_plan),
                 ('segment', args.segment),
                 ('observation', args.observation),
                 ('visit', args.visit)):
        if v is not None:
            globals_[k] = v

    def _merge(spec_str):
        spec = parse_spec(spec_str)
        # Per-layer values override the globals.
        return {**globals_, **spec}

    specs = {'blue': _merge(args.blue),
             'green': _merge(args.green),
             'red': _merge(args.red)}

    if args.mosaic:
        out_dir = (os.path.abspath(args.out_dir) if args.out_dir
                   else _default_out_dir(_spec_dir_name(specs) + '_mosaic'))
        os.makedirs(out_dir, exist_ok=True)
        print(f'[rgb] output dir: {out_dir}', file=sys.stderr)
        run_mosaic(
            specs, out_dir, workers=args.workers,
            from_cache=args.from_cache,
            save_png=_resolve_png_path(args.save_png, out_dir,
                                       'rgb_mosaic.png'),
            save_fits=args.save_fits,
            no_ds9=args.no_ds9,
            png_wcs=args.png_wcs,
        )
        return

    if args.save_fits or args.png_wcs:
        ap.error('--save-fits and --png-wcs only apply with --mosaic')

    sca = args.sca
    out_dir = (os.path.abspath(args.out_dir) if args.out_dir
               else _default_out_dir(_spec_dir_name(specs, sca=sca)))
    os.makedirs(out_dir, exist_ok=True)
    print(f'[rgb] output dir: {out_dir}', file=sys.stderr)
    png_path = _resolve_png_path(args.save_png, out_dir,
                                 f'rgb_sca{sca:02d}.png')

    layers = {}

    if args.from_cache:
        # For cache mode we need channel filenames — reconstruct from spec.
        # We don't know the actual filter name without streaming, so also
        # accept `filter=Fxxx` in the spec (optional).
        for channel, spec in specs.items():
            obs = int(spec.get('observation', 0))
            expn = int(spec['exposure'])
            filt = str(spec.get('optical_element',
                                spec.get('filter', 'UNKNOWN'))).upper()
            candidate = os.path.join(
                out_dir,
                cache_name(channel, {'obs': obs, 'exp': expn, 'filter': filt},
                           sca),
            )
            if not os.path.exists(candidate):
                # Fall back to glob: any file matching channel + sca.
                import glob
                pattern = os.path.join(
                    out_dir, f'{channel}_*_sca{sca:02d}.fits',
                )
                matches = sorted(glob.glob(pattern))
                if not matches:
                    raise FileNotFoundError(
                        f'--from-cache: no {channel} FITS matches {pattern}. '
                        f'Run without --from-cache to generate them.'
                    )
                candidate = matches[0]
            print(f'[rgb] loading {candidate}', file=sys.stderr)
            with fits.open(candidate) as hdul:
                data = np.asarray(hdul[0].data, dtype=np.float32)
                hdr = hdul[0].header
            layers[channel] = {'aligned': data, 'hdr': hdr,
                               'obs': obs, 'exp': expn, 'filter': filt}
        blue = layers['blue']
        if png_path:
            _sca_png(layers, png_path)
        if not args.no_ds9:
            _push_ds9(_run_ds9, layers, blue, file_outputs=bool(png_path))
        print('[rgb] done (cache).', file=sys.stderr)
        return

    # Stream all three layers.
    for channel, spec in specs.items():
        data, hdr, filt, obs, expn, visit_id = fetch_sca(spec, sca)
        layers[channel] = {'data': data, 'hdr': hdr, 'filter': filt,
                           'obs': obs, 'exp': expn, 'visit_id': visit_id}

    # Pure pixel-space alignment. Roman commissioning WCS isn't trustworthy
    # yet (the whole point of this stack is helping calibrate it), so we
    # DON'T reproject — every SCA n comes off the same detector at the same
    # 4088×4088 pixel grid, and the small dither between exposures shows up
    # as an integer-ish (dy, dx) shift we can measure directly.
    blue = layers['blue']
    blue['aligned'] = blue['data']

    for channel in ('green', 'red'):
        lyr = layers[channel]
        lyr['aligned'] = align_to_ref(lyr['data'], blue['aligned'],
                                      label=channel)

    # Write aligned FITS (all share the blue WCS header).
    for channel in ('red', 'green', 'blue'):
        lyr = layers[channel]
        path = os.path.join(out_dir, cache_name(channel, lyr, sca))
        write_fits(path, lyr['aligned'], blue['hdr'])

    # PNG before DS9, so a missing DS9 can't cost us the file.
    if png_path:
        _sca_png(layers, png_path)
    pushed = False
    if not args.no_ds9:
        pushed = _push_ds9(_run_ds9, layers, blue,
                           file_outputs=bool(png_path))
    print(f'[rgb] done. FITS written to {out_dir}'
          f'{". RGB loaded in DS9" if pushed else ""}', file=sys.stderr)


if __name__ == '__main__':
    main()
