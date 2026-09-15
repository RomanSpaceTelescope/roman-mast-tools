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


def fetch_sca(spec: dict, sca: int):
    """Return (data, wcs_header, filter_name, obs, expn, visit_id)."""
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

    dm_dict = stream_materialized(exp, res.missions, scas=[sca], max_workers=1)
    dm = dm_dict.get(sca)
    if dm is None:
        raise RuntimeError(f'SCA {sca} missing for {filters} '
                           f'exp={exposure_number}')
    try:
        data = np.asarray(dm.data[...], dtype=np.float32)
        gwcs = dm.meta.wcs
        wcs_hdr = gwcs.to_fits_sip(bounding_box=gwcs.bounding_box, degree=4)
        filt = str(exp.optical_element)
    finally:
        close_streams(dm_dict)
    obs = int(getattr(exp, 'observation', 0) or 0)
    return data, wcs_hdr, filt, obs, exposure_number, exp.visit_id


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


def _clean_for_xcorr(img):
    arr = np.array(img, dtype=np.float32)
    finite = np.isfinite(arr)
    if finite.any():
        _, med, _ = sigma_clipped_stats(arr[finite])
    else:
        med = 0.0
    return np.where(finite, arr - med, 0.0)


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
    d.set('zoom to fit')


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--sca', type=int, required=True, help='SCA number (1-18)')
    ap.add_argument('--program', type=int, default=None,
                    help='Global program ID inherited by all three layers '
                         '(override per-layer via the spec).')
    ap.add_argument('--pass', dest='pass_', type=int, default=None,
                    help='Global pass number inherited by all three layers '
                         '(override per-layer via the spec).')
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
    args = ap.parse_args()

    sca = args.sca
    out_dir = os.path.abspath(args.out_dir or f'rgb_sca{sca:02d}')
    os.makedirs(out_dir, exist_ok=True)

    globals_ = {}
    if args.program is not None:
        globals_['program'] = args.program
    if args.pass_ is not None:
        globals_['pass_'] = args.pass_

    def _merge(spec_str):
        spec = parse_spec(spec_str)
        # Per-layer values override the globals.
        return {**globals_, **spec}

    specs = {'blue': _merge(args.blue),
             'green': _merge(args.green),
             'red': _merge(args.red)}

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

    # Bounded-search cross-correlation. Real dithers between these visits are
    # sub-arcminute (small tens of px), so we compute the full FFT-based
    # cross-correlation and look for the peak within a ±SEARCH_RADIUS window
    # around zero. This avoids unbounded phase_cross_correlation locking
    # onto spurious peaks 1000+ px away when the sky content differs
    # strongly across filters.
    SEARCH_RADIUS_PX = 10

    def _bounded_shift(ref_full, mov_full, search_radius=SEARCH_RADIUS_PX):
        ref = _clean_for_xcorr(ref_full)
        mov = _clean_for_xcorr(mov_full)
        # Normalized cross-correlation via FFT.
        corr = fftconvolve(ref, mov[::-1, ::-1], mode='same')
        cy, cx = corr.shape[0] // 2, corr.shape[1] // 2
        # Restrict to the ±search_radius window around zero shift.
        y0, y1 = cy - search_radius, cy + search_radius + 1
        x0, x1 = cx - search_radius, cx + search_radius + 1
        window = corr[y0:y1, x0:x1]
        py, px = np.unravel_index(np.argmax(window), window.shape)
        dy = (py + y0) - cy
        dx = (px + x0) - cx
        return float(dy), float(dx), float(window.max())

    for channel in ('green', 'red'):
        lyr = layers[channel]
        raw = lyr['data']
        finite_in = np.isfinite(raw)
        print(f'[rgb] {channel} raw: finite={finite_in.sum()}/{raw.size} '
              f'({100.0*finite_in.mean():.1f}%)', file=sys.stderr)
        dy, dx, peak = _bounded_shift(blue['aligned'], raw)
        print(f'[rgb] {channel} shift (dy, dx) = ({dy:+.1f}, {dx:+.1f}) px  '
              f'peak={peak:.3g}  (searched ±{SEARCH_RADIUS_PX} px)',
              file=sys.stderr)
        # ndi_shift order>=1 with NaN input propagates NaN across a wide
        # region via spline interpolation. Zero-fill the NaN pixels first,
        # then apply the shift, then re-NaN the outside-frame footprint.
        clean = np.where(finite_in, raw, 0.0)
        shifted = ndi_shift(
            clean, shift=(dy, dx), order=3,
            mode='constant', cval=np.nan,
        )
        # Also shift the finite mask (nearest-neighbour) so we can restore
        # NaN in what came from outside the source frame.
        finite_shifted = ndi_shift(
            finite_in.astype(np.float32), shift=(dy, dx), order=0,
            mode='constant', cval=0.0,
        ) > 0.5
        shifted[~finite_shifted] = np.nan
        lyr['aligned'] = shifted.astype(np.float32)
        print(f'[rgb] {channel} aligned: finite='
              f'{int(np.isfinite(shifted).sum())}/{shifted.size} '
              f'({100.0*np.isfinite(shifted).mean():.1f}%)',
              file=sys.stderr)

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
