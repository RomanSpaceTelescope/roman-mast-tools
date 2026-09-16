"""Build an RGB composite of one SCA in DS9.

Each RGB channel is specified as a comma-separated key=value list of
`roman_mast.list_data` filters, plus `exposure=N` to pick which exposure
inside that query. Any filter accepted by roman_mast (program, pass,
execution_plan, segment, observation, visit, detector, optical_element,
data_level) can be used.

Example (program 1047 pass 1, three filters):

    python rgb_sca9_1047_p1.py --sca 7 \\
        --blue  program=1047,pass=1,observation=12,exposure=2 \\
        --green program=1047,pass=1,observation=5,exposure=2  \\
        --red   program=1047,pass=1,observation=9,exposure=2
"""
from __future__ import annotations

import argparse
import io
import os
import sys

import numpy as np
import pyds9
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import shift as ndi_shift
from scipy.signal import fftconvolve

from roman_mast import list_data, close_streams
from roman_fits import stream_materialized

DS9_TARGET = None  # None → pyds9 default; else the XPA target name

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
    hdu.writeto(path, overwrite=True)
    print(f'[rgb] wrote {path}', file=sys.stderr)


def _fits_bytes(data: np.ndarray, hdr: fits.Header) -> bytes:
    hdulist = fits.HDUList([fits.PrimaryHDU(data=data.astype(np.float32),
                                            header=hdr)])
    buf = io.BytesIO()
    hdulist.writeto(buf)
    return buf.getvalue()


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


def _run_ds9(layers, blue):
    """Push the three aligned layers into DS9 as an RGB frame."""
    GAIN_SIGMA = {'blue': 30.0, 'green': 25.0, 'red': 20.0}
    LOW_SIGMA = 1.0

    print('[rgb] connecting to DS9 via pyds9', file=sys.stderr)
    d = pyds9.DS9(target=DS9_TARGET) if DS9_TARGET else pyds9.DS9()
    d.set('frame delete all')
    d.set('frame new rgb')
    for channel in ('red', 'green', 'blue'):
        lyr = layers[channel]
        d.set(f'rgb channel {channel}')
        d.set('fits', _fits_bytes(lyr['aligned'], blue['hdr']))

        arr = lyr['aligned']
        finite = np.isfinite(arr)
        if finite.any():
            _, med, std = sigma_clipped_stats(arr[finite])
        else:
            med, std = 0.0, 1.0
        lo = float(med - LOW_SIGMA * std)
        hi = float(med + GAIN_SIGMA[channel] * std)
        print(f'[rgb] {channel}: bkg={med:.4g}  rms={std:.4g}  '
              f'limits=[{lo:.4g}, {hi:.4g}]  asinh', file=sys.stderr)
        d.set(f'scale limits {lo} {hi}')
        d.set('scale asinh')
    d.set('rgb channel red')
    d.set('wcs align no')  # image-aligned display; WCS not yet trustworthy
    d.set('zoom to fit')


ALL_SCAS = list(range(1, 19))


def _stretch_per_channel(all_arrays_by_channel):
    """Compute one (lo, hi) per channel from pooled sigma-clipped stats.

    all_arrays_by_channel : {channel: [aligned_ndarray, ...]}
        Concatenated across SCAs, so the mosaic gets a uniform stretch.
    """
    GAIN_SIGMA = {'blue': 30.0, 'green': 25.0, 'red': 20.0}
    LOW_SIGMA = 1.0
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
    buf = io.BytesIO()
    hdul.writeto(buf)
    data = buf.getvalue()
    if out_path is not None:
        with open(out_path, 'wb') as f:
            f.write(data)
        print(f'[rgb] wrote {out_path} ({len(data)/1e6:.0f} MB, '
              f'{len(hdul)-1} SCAs)', file=sys.stderr)
    return data


def _headless_mosaic_png(sca_layers, limits, ref_hdrs, png_path: str):
    """Render an RGB PNG of the focal-plane mosaic without DS9.

    Uses reproject.mosaicking to reproject the 18 per-SCA aligned arrays
    of each channel onto a single common celestial WCS, applies the same
    asinh + per-channel limits used in DS9, stacks R/G/B, and saves via
    matplotlib.
    """
    from reproject import reproject_interp
    from reproject.mosaicking import (
        find_optimal_celestial_wcs, reproject_and_coadd,
    )
    import matplotlib.pyplot as plt

    def _asinh_scale(arr, lo, hi):
        # Match DS9's asinh: normalize to [0, 1] then apply asinh stretch.
        x = (arr - lo) / max(hi - lo, 1e-12)
        x = np.clip(x, 0.0, 1.0)
        return np.arcsinh(10.0 * x) / np.arcsinh(10.0)

    # Common WCS built from blue's per-SCA WCS headers.
    inputs_blue = [(sca_layers[sca]['blue'], WCS(ref_hdrs[sca]))
                   for sca in sorted(sca_layers)
                   if 'blue' in sca_layers[sca]]
    common_wcs, common_shape = find_optimal_celestial_wcs(inputs_blue)
    print(f'[rgb] headless mosaic grid: shape={common_shape} '
          f'CRVAL={common_wcs.wcs.crval}', file=sys.stderr)

    rgb = np.zeros((common_shape[0], common_shape[1], 3), dtype=np.float32)
    for idx, channel in enumerate(('red', 'green', 'blue')):
        inputs = [(sca_layers[sca][channel], WCS(ref_hdrs[sca]))
                  for sca in sorted(sca_layers)
                  if channel in sca_layers[sca]]
        mosaic, _ = reproject_and_coadd(
            inputs, common_wcs, shape_out=common_shape,
            reproject_function=reproject_interp,
            combine_function='mean',
        )
        lo, hi = limits[channel]
        rgb[:, :, idx] = _asinh_scale(mosaic, lo, hi)

    # Matplotlib expects sky-north-up. Astropy WCS puts DEC increasing in
    # array Y; flip vertically so the PNG looks like DS9's default view.
    rgb = np.flipud(rgb)
    fig, ax = plt.subplots(figsize=(12, 12), dpi=150)
    ax.imshow(rgb, origin='upper')
    ax.set_axis_off()
    fig.savefig(png_path, bbox_inches='tight', pad_inches=0,
                facecolor='black')
    plt.close(fig)
    print(f'[rgb] wrote headless PNG {png_path}', file=sys.stderr)


def _run_ds9_mosaic(sca_layers, limits, ref_hdrs, out_dir=None,
                    save_png=None):
    """Push a single WCS-stitched RGB mosaic into DS9.

    Builds three MEFs (one per channel, each holding all 18 SCAs as
    extensions with SIP WCS headers) and pipes each into DS9's
    `fits mosaicimage wcs` command on the matching rgb channel of one
    RGB frame. Every channel of the frame becomes a full-focal-plane
    mosaic, and DS9 combines the three into a single color view.

    Optionally also writes the three MEFs to `out_dir/mosaic_{rgb}.fits`
    so downstream tools can inspect them.
    """
    print(f'[rgb] connecting to DS9 via pyds9', file=sys.stderr)
    d = pyds9.DS9(target=DS9_TARGET) if DS9_TARGET else pyds9.DS9()
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

    if save_png:
        # Try a few DS9 save variants — behaviour depends on DS9 version
        # and whether the DS9 window is realized (some builds fail
        # `saveimage png` on headless / offscreen setups).
        import time
        abs_path = os.path.abspath(save_png)
        os.makedirs(os.path.dirname(abs_path) or '.', exist_ok=True)
        # Give DS9 a moment to finish rendering before we ask for the image.
        d.set('update now')
        time.sleep(0.5)
        variants = [
            f'saveimage png {abs_path}',
            f'saveimage {abs_path} png',
            f'export png {abs_path}',
            f'export {abs_path} png',
        ]
        saved = False
        for cmd in variants:
            try:
                d.set(cmd)
                if os.path.exists(abs_path) and os.path.getsize(abs_path) > 0:
                    print(f'[rgb] wrote DS9 PNG via `{cmd}` → {abs_path}',
                          file=sys.stderr)
                    saved = True
                    break
                else:
                    print(f'[rgb] DS9 `{cmd}` returned OK but no file at '
                          f'{abs_path}', file=sys.stderr)
            except Exception as e:
                print(f'[rgb] DS9 `{cmd}` failed: '
                      f'{type(e).__name__}: {e}', file=sys.stderr)
        if not saved:
            print('[rgb] WARNING: all DS9 saveimage/export variants failed. '
                  'Use --headless-png for a matplotlib render instead.',
                  file=sys.stderr)


def run_mosaic(specs, out_dir, *, workers=8, from_cache=False,
               save_png=None, headless_png=None):
    """Stream all 18 SCAs for each of the three channels, align per-SCA, push
    to DS9 as one RGB frame per SCA."""
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
        _run_ds9_mosaic(sca_layers, limits, ref_hdrs, out_dir=out_dir,
                        save_png=save_png)
        if headless_png:
            _headless_mosaic_png(sca_layers, limits, ref_hdrs, headless_png)
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
    _run_ds9_mosaic(sca_layers, limits, ref_hdrs, out_dir=out_dir,
                    save_png=save_png)
    if headless_png:
        _headless_mosaic_png(sca_layers, limits, ref_hdrs, headless_png)
    print(f'[rgb] mosaic done. FITS in {out_dir}', file=sys.stderr)


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
                         'independently (blue as pivot), and loads 18 RGB '
                         'frames into DS9 locked by WCS.')
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
                    help='Output directory for aligned FITS '
                         '(default: rgb_sca<N>)')
    ap.add_argument('--from-cache', action='store_true',
                    help='Skip MAST streaming + alignment; load the aligned '
                         'FITS files from --out-dir and push to DS9.')
    ap.add_argument('--save-png', default=None,
                    help='After loading in DS9, export the current view to '
                         'PNG via DS9 (WYSIWYG of the DS9 display). Path is '
                         'either absolute or relative to --out-dir.')
    ap.add_argument('--headless-png', default=None,
                    help='Render a PNG without DS9 via matplotlib. Uses the '
                         'same asinh + per-channel limits as the DS9 push. '
                         'Path is absolute or relative to --out-dir.')
    args = ap.parse_args()

    def _resolve_png_path(p, base):
        if p is None:
            return None
        return p if os.path.isabs(p) else os.path.join(base, p)

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
            save_png=_resolve_png_path(args.save_png, out_dir),
            headless_png=_resolve_png_path(args.headless_png, out_dir),
        )
        return

    sca = args.sca
    out_dir = (os.path.abspath(args.out_dir) if args.out_dir
               else _default_out_dir(_spec_dir_name(specs, sca=sca)))
    os.makedirs(out_dir, exist_ok=True)
    print(f'[rgb] output dir: {out_dir}', file=sys.stderr)

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
        _run_ds9(layers, blue)
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

    _run_ds9(layers, blue)
    print(f'[rgb] done. RGB loaded in DS9. FITS written to {out_dir}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
