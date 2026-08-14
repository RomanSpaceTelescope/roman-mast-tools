#!/usr/bin/env python
"""Run aperture photometry across all SCAs of a Roman WFI exposure.

Reads a plain-text file of S3 URIs (one per line), streams each SCA in
parallel, and writes a combined per-source CSV and a per-SCA summary CSV.

Usage
-----
    # Build a URI file (one S3 path per line, blank lines / # comments OK):
    #   s3://stpubdata/roman/.../r0003201001001001004_0001_wfi01_f106_cal.asdf
    #   s3://stpubdata/roman/.../r0003201001001001004_0001_wfi02_f106_cal.asdf
    #   ...

    python roman_phot.py --uri-file my_exposure.txt

    # All outputs land in {visit_id}_{exposure_num}/ e.g.:
    #   r0003201001001001004_0001/sources.csv
    #   r0003201001001001004_0001/summary.csv
    #   r0003201001001001004_0001/histograms.png
    #   r0003201001001001004_0001/bkg_mosaic.png       (with --bkg-mosaic)
    #   r0003201001001001004_0001/source_mosaic.png    (with --bkg-mosaic)
    #   r0003201001001001004_0001/mosaic_data.npz      (with --bkg-mosaic)
    #   r0003201001001001004_0001/sca{NN}.csv          (with --per-sca)

    # Photometry tuning:
    python roman_phot.py --uri-file my_exposure.txt \\
        --fwhm 1.5 --det-sigma 10 --snr-threshold 5

    # Per-SCA CSVs as well:
    python roman_phot.py --uri-file my_exposure.txt --per-sca

    # Background + source-density mosaics:
    python roman_phot.py --uri-file my_exposure.txt --bkg-mosaic

    # Re-render mosaics from a previous run (no photometry re-run):
    python roman_phot.py --remake-mosaics r0003201001001001004_0001/mosaic_data.npz
"""

import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import vstack, Table
from astropy.stats import sigma_clipped_stats
from astropy.modeling import models, fitting
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.visualization import simple_norm
from photutils.segmentation import detect_sources
from scipy.ndimage import binary_dilation, generate_binary_structure, iterate_structure

from roman_view_sca import stream_sca, parse_filename, display_in_mpl, display_residuals_mpl
from photometry import run_aperture_photometry

# ---------------------------------------------------------------------------
# WFI focal-plane constants
# ---------------------------------------------------------------------------

# Roman WFI detector constants
_ROMAN_PIXEL_SCALE_MM = 0.01      # mm per pixel (focal-plane scale)
_ROMAN_SCA_FULL_SIZE  = 4096      # full detector size before reference-pixel removal
_ROMAN_REF_PIX        = 4         # reference pixels removed from each edge

# Focal-plane layout: SCA number → (x_center_mm, y_center_mm, rotation_deg)
# Coordinates are focal-plane mm; rotation 180° means the SCA is flipped
# relative to the focal-plane axes.  Source: MPA_SCA_info ECSV (width=40.88 mm).
_WFI_SCA_LAYOUT = {
     1: (  -22.14,   12.15, 180.0),
     2: (  -22.29,  -37.03, 180.0),
     3: (  -22.44,  -82.06,   0.0),
     4: (  -66.42,   20.90, 180.0),
     5: (  -66.92,  -28.28, 180.0),
     6: (  -67.42,  -73.06,   0.0),
     7: ( -110.70,   42.20, 180.0),
     8: ( -111.48,   -6.98, 180.0),
     9: ( -112.64,  -51.06,   0.0),
    10: (   22.14,   12.15, 180.0),
    11: (   22.29,  -37.03, 180.0),
    12: (   22.44,  -82.06,   0.0),
    13: (   66.42,   20.90, 180.0),
    14: (   66.92,  -28.28, 180.0),
    15: (   67.42,  -73.06,   0.0),
    16: (  110.70,   42.20, 180.0),
    17: (  111.48,   -6.98, 180.0),
    18: (  112.64,  -51.06,   0.0),
}


# ---------------------------------------------------------------------------
# Per-SCA photometry
# ---------------------------------------------------------------------------

def phot_one_sca(uri, filename, *, phot_kwargs, bkg_kwargs=None, png_path=None,
                 image_mosaic=False):
    """Stream one SCA and run aperture photometry.

    Returns a dict with keys: sca, detector, table, stats, sp, bkg_map, thumb.
    bkg_map is None when bkg_kwargs is None; thumb is None when image_mosaic is False.
    Returns None on failure (logs a warning).
    """
    try:
        data, detector, _wcs_hdr, dq = stream_sca(uri, filename)
    except Exception as exc:
        print(f'[roman_phot] WARNING: failed to stream {filename}: {exc}', file=sys.stderr)
        return None

    parsed = parse_filename(filename)
    sca_num = parsed['sca_num']

    try:
        sources, phot_table, sp, stats = run_aperture_photometry(
            data, dq=dq, **phot_kwargs
        )
    except Exception as exc:
        print(f'[roman_phot] WARNING: photometry failed for SCA {sca_num}: {exc}', file=sys.stderr)
        return None

    if phot_table is not None and len(phot_table) > 0:
        phot_table['sca'] = sca_num
        phot_table['detector'] = detector
        # Move sca/detector to the front
        col_order = ['sca', 'detector'] + [c for c in phot_table.colnames
                                            if c not in ('sca', 'detector')]
        phot_table = phot_table[col_order]

    snr_arr = phot_table['snr'] if (phot_table is not None and
                                     len(phot_table) > 0 and
                                     'snr' in phot_table.colnames) else None

    print(
        f'[roman_phot] SCA {sca_num:02d} ({detector}): '
        f'bkg={stats["bkg_level"]:.3g}  rms={stats["bkg_rms_median"]:.3g}  '
        f'n_sources={stats["n_sources"]}' +
        (f'  snr=[{float(snr_arr.min()):.1f}, {float(np.median(snr_arr)):.1f}, {float(snr_arr.max()):.1f}]'
         if snr_arr is not None and len(snr_arr) > 0 else ''),
        file=sys.stderr,
    )

    bkg_map = None
    if bkg_kwargs is not None:
        try:
            bkg_map = background_map_one_sca(data, dq, data_sub=stats['data_sub'], **bkg_kwargs)
        except Exception as exc:
            print(f'[roman_phot] WARNING: background map failed for SCA {sca_num}: {exc}',
                  file=sys.stderr)

    if png_path is not None:
        try:
            display_in_mpl(
                data, detector, dq=dq,
                title=f'{detector}  {filename}',
                sources=sources, phot_table=phot_table, sp=sp, stats=stats,
                save_path=png_path,
            )
            print(f'[roman_phot] SCA {sca_num:02d} image -> {png_path}', file=sys.stderr)
        except Exception as exc:
            print(f'[roman_phot] WARNING: image PNG failed for SCA {sca_num}: {exc}',
                  file=sys.stderr)

        bkg_fit = stats.get('bkg_fit')
        data_sub = stats.get('data_sub')
        if data_sub is not None:
            resid_path = png_path.replace('.png', '_bkg.png')
            try:
                display_residuals_mpl(
                    data_sub, detector,
                    bkg_fit=bkg_fit,
                    bkg_level=stats['bkg_level'],
                    residual_rms=stats['bkg_rms_median'],
                    poly_degree=stats.get('poly_degree'),
                    dq=dq,
                    save_path=resid_path,
                )
                print(f'[roman_phot] SCA {sca_num:02d} bkg/residuals -> {resid_path}',
                      file=sys.stderr)
            except Exception as exc:
                print(f'[roman_phot] WARNING: bkg PNG failed for SCA {sca_num}: {exc}',
                      file=sys.stderr)

    thumb = data[::2, ::2].astype(np.float32) if image_mosaic else None

    return {
        'sca': sca_num,
        'detector': detector,
        'table': phot_table,
        'stats': stats,
        'sp': sp,
        'bkg_map': bkg_map,
        'thumb': thumb,
    }


# ---------------------------------------------------------------------------
# Background mapping
# ---------------------------------------------------------------------------

def background_map_one_sca(data, dq, *, superpixel=512, mask_sigma=1.5,
                            dilate_radius=20, bkg_poly_degree=3, data_sub=None):
    """Compute a superpixel background map for one SCA.

    Steps:
      1. Fit and subtract a 2D polynomial background (skipped if data_sub provided).
      2. Mask all sources (point and extended) via a segmentation map with
         binary dilation so halos are fully covered.
      3. Replace masked pixels with NaN.
      4. Bin into superpixels using nanmedian.

    Returns a 2D float array of shape (ny//superpixel, nx//superpixel), or
    None if the image is too small.
    """
    ny, nx = data.shape
    mask = (dq != 0) if dq is not None else np.zeros_like(data, dtype=bool)

    if data_sub is None:
        # Polynomial background fit (same logic as run_aperture_photometry)
        poly_init = models.Polynomial2D(degree=bkg_poly_degree)
        fitter = fitting.LinearLSQFitter()
        ds = 4
        y_ds, x_ds = np.mgrid[0:ny:ds, 0:nx:ds]
        valid = ~mask[::ds, ::ds]
        if np.any(valid):
            data_ds = data[::ds, ::ds]
            poly = fitter(poly_init, x_ds[valid], y_ds[valid], data_ds[valid])
            y_full, x_full = np.mgrid[:ny, :nx]
            bkg_fit = poly(x_full, y_full)
        else:
            bkg_fit = np.zeros_like(data, dtype=float)
        data_sub = data - bkg_fit

    # RMS of background residual (sigma-clipped, excluding bad pixels)
    residuals = data_sub[~mask] if np.any(~mask) else data_sub.ravel()
    residuals = residuals[np.isfinite(residuals)]
    _, _, rms = sigma_clipped_stats(residuals, sigma=5.0)

    # Build source mask via segmentation
    src_mask = np.zeros_like(mask)
    try:
        kernel = Gaussian2DKernel(x_stddev=2.0)
        conv = convolve(data_sub, kernel, mask=mask, nan_treatment='fill', preserve_nan=False)
        segmap = detect_sources(conv, mask_sigma * rms, n_pixels=5)
        if segmap is not None:
            # Single-pass dilation with a pre-expanded disk avoids dilate_radius
            # serial GIL-holding iterations over the full array.
            disk = iterate_structure(generate_binary_structure(2, 1), dilate_radius)
            src_mask = binary_dilation(segmap.data > 0, structure=disk)
    except Exception as exc:
        print(f'[roman_phot] WARNING: source masking skipped: {exc}', file=sys.stderr)

    # Apply combined mask as NaN
    residual = data_sub.astype(float)
    residual[mask | src_mask] = np.nan

    # Bin into superpixels aligned to the original 4096×4096 detector grid.
    # Reference pixels shift all data coordinates by -_ROMAN_REF_PIX; edge
    # superpixels end up with _ROMAN_REF_PIX fewer rows/columns, which is fine.
    n_chunks = _ROMAN_SCA_FULL_SIZE // superpixel  # 8 for superpixel=512

    orig_bounds = [i * superpixel for i in range(n_chunks + 1)]
    row_bounds = [max(0, min(ny, b - _ROMAN_REF_PIX)) for b in orig_bounds]
    col_bounds = [max(0, min(nx, b - _ROMAN_REF_PIX)) for b in orig_bounds]

    binned = np.full((n_chunks, n_chunks), np.nan)
    for i in range(n_chunks):
        r0, r1 = row_bounds[i], row_bounds[i + 1]
        if r1 <= r0:
            continue
        for j in range(n_chunks):
            c0, c1 = col_bounds[j], col_bounds[j + 1]
            if c1 <= c0:
                continue
            chunk = residual[r0:r1, c0:c1]
            binned[i, j] = np.nanmedian(chunk)

    n_src = int(np.sum(src_mask))
    n_tot = src_mask.size
    print(
        f'[roman_phot] bkg map: {n_chunks}×{n_chunks} superpixels  '
        f'source mask={100*n_src/n_tot:.1f}%  rms={rms:.4g}',
        file=sys.stderr,
    )
    return binned



def make_image_mosaic_png(sca_thumbs, out_path, *, title=None):
    """Render a WFI focal-plane image mosaic from decimated SCA thumbnails.

    Parameters
    ----------
    sca_thumbs : dict
        SCA number → 2D float32 array decimated by 4 in each axis (or None).
    out_path : str
        Destination PNG path.
    title : str or None
        Figure title.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    thumbs = {k: v for k, v in sca_thumbs.items() if v is not None}
    if not thumbs:
        print('[roman_phot] WARNING: no image thumbnails to render; skipping image mosaic',
              file=sys.stderr)
        return

    all_data = np.concatenate([v.ravel() for v in thumbs.values()])
    norm = simple_norm(all_data, 'asinh', vmin=0.5, vmax=4)

    half = _ROMAN_SCA_FULL_SIZE * _ROMAN_PIXEL_SCALE_MM / 2  # 20.48 mm

    all_cx = [cx for cx, _, _ in _WFI_SCA_LAYOUT.values()]
    all_cy = [cy for _, cy, _ in _WFI_SCA_LAYOUT.values()]
    pad = 5.0
    x_lo = min(all_cx) - half - pad
    x_hi = max(all_cx) + half + pad
    y_lo = min(all_cy) - half - pad
    y_hi = max(all_cy) + half + pad

    aspect = (x_hi - x_lo) / (y_hi - y_lo)
    fig_w = 14.0
    fig_h = fig_w / aspect
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor='#1a1a1a')
    ax.set_facecolor('#1a1a1a')
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect('equal')

    for sca_num, (cx_mm, cy_mm, rot) in _WFI_SCA_LAYOUT.items():
        tile = thumbs.get(sca_num)
        has_data = tile is not None

        x0 = cx_mm - half
        x1 = cx_mm + half
        y0 = cy_mm - half
        y1 = cy_mm + half

        if has_data:
            if rot == 180.0:
                tile = np.rot90(tile, 2)
            ax.imshow(tile, extent=[x0, x1, y0, y1],
                      origin='upper', cmap='gray', norm=norm,
                      interpolation='nearest', aspect='auto')

        rect = Rectangle(
            (x0, y0), x1 - x0, y1 - y0,
            linewidth=0.6,
            edgecolor='#555555' if has_data else '#333333',
            facecolor='none',
        )
        ax.add_patch(rect)
        if not has_data:
            ax.text(cx_mm, cy_mm, f'{sca_num:02d}',
                    ha='center', va='center', fontsize=7,
                    color='#555555', fontweight='bold')

    ax.set_title(title or 'Roman WFI — image mosaic (::2)',
                 color='white', fontsize=12, pad=10)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'[roman_phot] image mosaic -> {out_path}', file=sys.stderr)


def make_bkg_mosaic_png(sca_maps, out_path, *, superpixel=512, title=None,
                        pct_lo=2, pct_hi=98):
    """Render a WFI focal-plane background mosaic and save to a PNG.

    Parameters
    ----------
    sca_maps : dict
        Mapping of SCA number (int) to 2D binned background array (or None).
    out_path : str
        Destination PNG path.
    superpixel : int
        Superpixel size used when creating the maps (for scale conversion).
    title : str or None
        Figure title.
    pct_lo, pct_hi : float
        Percentile bounds for the symmetric colour stretch.
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    from matplotlib.patches import Rectangle

    # Physical size of one superpixel in focal-plane mm
    sp_mm = superpixel * _ROMAN_PIXEL_SCALE_MM   # e.g. 512 × 0.01 = 5.12 mm

    # Infer tile shape from available maps
    tile_shapes = [v.shape for v in sca_maps.values() if v is not None]
    if not tile_shapes:
        print(f'[roman_phot] WARNING: no background maps to render; skipping PNG', file=sys.stderr)
        return
    tile_h, tile_w = tile_shapes[0]

    # Physical half-size of one SCA tile in mm
    half_w_mm = tile_w * sp_mm / 2
    half_h_mm = tile_h * sp_mm / 2

    # Axes limits in focal-plane mm (y+ = up, matching sky orientation)
    all_cx = [cx for cx, _, _ in _WFI_SCA_LAYOUT.values()]
    all_cy = [cy for _, cy, _ in _WFI_SCA_LAYOUT.values()]
    pad = 5.0   # mm padding around outermost SCAs
    x_lo = min(all_cx) - half_w_mm - pad
    x_hi = max(all_cx) + half_w_mm + pad
    y_lo = min(all_cy) - half_h_mm - pad
    y_hi = max(all_cy) + half_h_mm + pad

    # Colour stretch: symmetric around zero (diverging), from all available tiles
    all_vals = np.concatenate([v.ravel() for v in sca_maps.values()
                                if v is not None and np.any(np.isfinite(v))])
    all_vals = all_vals[np.isfinite(all_vals)]
    if len(all_vals) == 0:
        print(f'[roman_phot] WARNING: background mosaic is entirely NaN; skipping PNG', file=sys.stderr)
        return
    abs_lim = max(abs(np.percentile(all_vals, pct_lo)), abs(np.percentile(all_vals, pct_hi)))
    norm = mcolors.Normalize(vmin=-abs_lim, vmax=abs_lim)

    aspect = (x_hi - x_lo) / (y_hi - y_lo)
    fig_w = 12.0
    fig_h = fig_w / aspect
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor='#1a1a1a')
    ax.set_facecolor('#1a1a1a')
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect('equal')

    # Draw each SCA tile as its own imshow in mm coordinates — no integer
    # rounding, so physical spacing between SCAs is preserved exactly.
    for sca_num, (cx_mm, cy_mm, rot) in _WFI_SCA_LAYOUT.items():
        tile = sca_maps.get(sca_num)
        has_data = tile is not None

        x0 = cx_mm - half_w_mm
        x1 = cx_mm + half_w_mm
        y0 = cy_mm - half_h_mm
        y1 = cy_mm + half_h_mm

        if has_data:
            if rot == 180.0:
                tile = np.rot90(tile, 2)
            # extent=[left, right, bottom, top]; origin='upper' maps row 0 to top
            ax.imshow(
                tile, extent=[x0, x1, y0, y1],
                origin='upper', cmap='RdBu_r', norm=norm,
                interpolation='nearest', aspect='auto',
            )

        rect = Rectangle(
            (x0, y0), x1 - x0, y1 - y0,
            linewidth=0.6,
            edgecolor='white' if has_data else '#555555',
            facecolor='none',
            alpha=0.8 if has_data else 0.4,
        )
        ax.add_patch(rect)
        ax.text(cx_mm, cy_mm, f'{sca_num:02d}',
                ha='center', va='center', fontsize=7,
                color='white' if has_data else '#777777',
                alpha=0.8 if has_data else 0.4,
                fontweight='bold')

    sm = cm.ScalarMappable(norm=norm, cmap='RdBu_r')
    cbar = fig.colorbar(sm, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label('Background residual (DN/s)', color='white', fontsize=10)
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')

    ax.set_title(title or 'Roman WFI — source-masked background mosaic',
                 color='white', fontsize=12, pad=10)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'[roman_phot] background mosaic -> {out_path}', file=sys.stderr)


def make_source_dot_mosaic_png(sources_csv_path, out_path, *, title=None):
    """Render all detected sources as dots in the WFI focal-plane layout.

    Parameters
    ----------
    sources_csv_path : str
        Path to the sources CSV produced by roman_phot (must have sca,
        x_centroid, y_centroid columns).
    out_path : str
        Destination PNG path.
    title : str or None
        Figure title.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    try:
        sources = Table.read(sources_csv_path, format='ascii.csv')
    except Exception as exc:
        print(f'[roman_phot] WARNING: cannot read sources for dot mosaic: {exc}', file=sys.stderr)
        return

    # Support both column naming conventions (photutils uses x_center/y_center;
    # detection catalog uses x_centroid/y_centroid).
    if 'x_center' in sources.colnames:
        sources.rename_column('x_center', 'x_centroid')
        sources.rename_column('y_center', 'y_centroid')
    required = {'sca', 'x_centroid', 'y_centroid'}
    if not required.issubset(sources.colnames):
        print(f'[roman_phot] WARNING: sources CSV missing columns {required - set(sources.colnames)}; '
              f'skipping dot mosaic', file=sys.stderr)
        return

    half = _ROMAN_SCA_FULL_SIZE * _ROMAN_PIXEL_SCALE_MM / 2  # 20.48 mm

    all_cx = [cx for cx, _, _ in _WFI_SCA_LAYOUT.values()]
    all_cy = [cy for _, cy, _ in _WFI_SCA_LAYOUT.values()]
    pad = 5.0
    x_lo = min(all_cx) - half - pad
    x_hi = max(all_cx) + half + pad
    y_lo = min(all_cy) - half - pad
    y_hi = max(all_cy) + half + pad

    aspect = (x_hi - x_lo) / (y_hi - y_lo)
    fig_w = 14.0
    fig_h = fig_w / aspect
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor='#1a1a1a')
    ax.set_facecolor('#1a1a1a')
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect('equal')

    for sca_num, (cx_mm, cy_mm, rot) in _WFI_SCA_LAYOUT.items():
        mask = np.asarray(sources['sca']) == sca_num
        has_data = np.any(mask)

        x0 = cx_mm - half
        x1 = cx_mm + half
        y0 = cy_mm - half
        y1 = cy_mm + half

        if has_data:
            x_pix = np.asarray(sources['x_centroid'][mask], dtype=float)
            y_pix = np.asarray(sources['y_centroid'][mask], dtype=float)
            # Convert science pixel coords to focal-plane mm.
            # origin='upper' convention: col 0 → x0, row 0 → y1 (top).
            offset_x = (x_pix + _ROMAN_REF_PIX + 0.5) * _ROMAN_PIXEL_SCALE_MM
            offset_y = (y_pix + _ROMAN_REF_PIX + 0.5) * _ROMAN_PIXEL_SCALE_MM
            sign = -1 if rot == 180 else 1
            x_fp = cx_mm + sign * (offset_x - half)
            y_fp = cy_mm - sign * (offset_y - half)
            ax.scatter(x_fp, y_fp, s=0.5, c='#ffdd88', alpha=0.6,
                       linewidths=0, rasterized=True)

        rect = Rectangle(
            (x0, y0), x1 - x0, y1 - y0,
            linewidth=0.6,
            edgecolor='white' if has_data else '#555555',
            facecolor='none',
            alpha=0.8 if has_data else 0.4,
        )
        ax.add_patch(rect)
        ax.text(cx_mm, cy_mm, f'{sca_num:02d}',
                ha='center', va='center', fontsize=7,
                color='white' if has_data else '#777777',
                alpha=0.8 if has_data else 0.4,
                fontweight='bold')

    ax.set_title(title or f'Roman WFI — sources  ({len(sources):,} total)',
                 color='white', fontsize=12, pad=10)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'[roman_phot] source dot mosaic -> {out_path}', file=sys.stderr)


# ---------------------------------------------------------------------------
# Parallel fan-out
# ---------------------------------------------------------------------------

def phot_exposure(uri_filename_pairs, *, phot_kwargs=None, bkg_kwargs=None,
                  out_dir=None, max_workers=8, image_mosaic=False):
    """Run photometry on a list of (uri, filename) pairs in parallel.

    Returns a list of result dicts (sorted by SCA number), with None entries
    removed.
    """


    if phot_kwargs is None:
        phot_kwargs = {}

    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {}
        for uri, fn in uri_filename_pairs:
            meta = parse_filename(fn)
            sca_num = meta['sca_num']
            png_path = os.path.join(out_dir, f'sca{sca_num:02d}.png') if out_dir else None
            futures[pool.submit(phot_one_sca, uri, fn,
                                phot_kwargs=phot_kwargs, bkg_kwargs=bkg_kwargs,
                                png_path=png_path, image_mosaic=image_mosaic)] = fn
        for fut in as_completed(futures):
            res = fut.result()
            if res is not None:
                results[res['sca']] = res

    return [results[k] for k in sorted(results)]


def stream_image_mosaic(uri_filename_pairs, *, max_workers=8):
    """Stream all SCAs and return a dict of {sca_num: decimated thumbnail}."""
    def _stream_one(uri, fn):
        try:
            data, _det, _wcs, _dq = stream_sca(uri, fn)
            sca_num = parse_filename(fn)['sca_num']
            return sca_num, data[::2, ::2].astype(np.float32)
        except Exception as exc:
            print(f'[roman_phot] WARNING: failed to stream {fn}: {exc}', file=sys.stderr)
            return None

    thumbs = {}
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(_stream_one, uri, fn) for uri, fn in uri_filename_pairs]
        for fut in as_completed(futures):
            res = fut.result()
            if res is not None:
                thumbs[res[0]] = res[1]
    return thumbs


# ---------------------------------------------------------------------------
# URI file parsing
# ---------------------------------------------------------------------------

def load_uri_file(path):
    """Read a text file of S3 URIs. Returns list of (base_uri, filename) tuples."""
    pairs = []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            idx = line.rfind('/')
            if idx == -1:
                print(f'[roman_phot] WARNING: cannot parse URI (no /): {line!r}', file=sys.stderr)
                continue
            base_uri = line[:idx]
            filename = line[idx + 1:]
            pairs.append((base_uri, filename))
    return pairs


# ---------------------------------------------------------------------------
# Summary table builder
# ---------------------------------------------------------------------------

def build_summary(results):
    """Build an astropy Table with one row per SCA."""
    rows = []
    for r in results:
        sca = r['sca']
        det = r['detector']
        stats = r['stats']
        sp = r['sp']
        tbl = r['table']

        snr_min = snr_max = snr_median = float('nan')
        if tbl is not None and len(tbl) > 0 and 'snr' in tbl.colnames:
            snr = tbl['snr']
            snr_min    = float(snr.min())
            snr_max    = float(snr.max())
            snr_median = float(np.median(snr))

        rows.append({
            'sca':               sca,
            'detector':          det,
            'bkg_level':         stats['bkg_level'],
            'bkg_rms':           stats['bkg_rms_median'],
            'detection_threshold': stats['threshold'],
            'n_sources':         stats['n_sources'],
            'snr_min':           snr_min,
            'snr_median':        snr_median,
            'snr_max':           snr_max,
            'aperture_radius_pix': sp.aperture_radius if sp else float('nan'),
            'annulus_inner_pix': sp.annulus_inner    if sp else float('nan'),
            'annulus_outer_pix': sp.annulus_outer    if sp else float('nan'),
        })

    return Table(rows=rows) if rows else Table()


# ---------------------------------------------------------------------------
# Mosaic data persistence
# ---------------------------------------------------------------------------

def save_mosaic_data(sca_maps, out_path):
    """Save background maps to a compressed .npz file.

    Parameters
    ----------
    sca_maps : dict
        SCA number → 2D background array (or None).
    out_path : str
        Destination .npz path.
    """
    arrays = {f'bkg_{sca_num:02d}': arr
              for sca_num, arr in sca_maps.items() if arr is not None}
    np.savez_compressed(out_path, **arrays)
    print(f'[roman_phot] mosaic data -> {out_path}', file=sys.stderr)


def load_mosaic_data(path):
    """Load background maps previously saved with save_mosaic_data.

    Returns
    -------
    sca_maps : dict
    """
    data = np.load(path)
    sca_maps = {}
    for key in data.files:
        if key.startswith('bkg_'):
            sca_maps[int(key[4:])] = data[key]
    print(f'[roman_phot] loaded mosaic data from {path} ({len(sca_maps)} bkg maps)', file=sys.stderr)
    return sca_maps


# ---------------------------------------------------------------------------
# Histogram generation
# ---------------------------------------------------------------------------

def reject_outliers(data, method='iqr', iqr_mult=1.5):
    """Reject outliers using IQR or sigma clipping."""
    if method == 'iqr':
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        lower = q1 - iqr_mult * iqr
        upper = q3 + iqr_mult * iqr
        mask = (data >= lower) & (data <= upper)
    else:
        raise ValueError(f'Unknown outlier method: {method}')

    n_rejected = len(data) - np.sum(mask)
    return data[mask], n_rejected


def make_histograms_from_csv(csv_path, output_path, columns=None, bins=30, outlier_method='iqr'):
    """Generate histograms from a photometry CSV file.

    Parameters
    ----------
    csv_path : str
        Path to input CSV
    output_path : str
        Path to save output image
    columns : list of str, optional
        Columns to histogram (default: aperture_sum_0, snr)
    bins : int
        Number of bins
    outlier_method : str
        Outlier rejection method ('iqr' or 'none')
    """
    if columns is None:
        columns = ['aperture_sum_0', 'snr']

    try:
        table = Table.read(csv_path, format='ascii.csv')
    except Exception as exc:
        print(f'[roman_phot] WARNING: failed to read {csv_path} for histograms: {exc}',
              file=sys.stderr)
        return

    # Validate columns exist
    missing = [c for c in columns if c not in table.colnames]
    if missing:
        print(f'[roman_phot] WARNING: histogram columns not found: {missing}',
              file=sys.stderr)
        return

    n_cols = len(columns)
    n_rows = (n_cols + 1) // 2
    fig, axes = plt.subplots(n_rows, 2, figsize=(12, 4 * n_rows))
    if n_cols == 1:
        axes = np.array([[axes, axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)

    axes_flat = axes.flatten()

    for idx, colname in enumerate(columns):
        ax = axes_flat[idx]
        data = table[colname]

        # Filter non-finite
        valid_data = data[np.isfinite(data)]
        if len(valid_data) == 0:
            ax.text(0.5, 0.5, f'No valid data in {colname}',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(colname)
            continue

        # Reject outliers
        plot_data = valid_data
        n_rejected = 0
        if outlier_method == 'iqr':
            plot_data, n_rejected = reject_outliers(valid_data, method='iqr')

        if len(plot_data) == 0:
            ax.text(0.5, 0.5, f'No data left after outlier rejection',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(colname)
            continue

        ax.hist(plot_data, bins=bins, edgecolor='black', alpha=0.7)
        ax.set_xlabel(colname)
        ax.set_ylabel('Count')
        title_str = f'{colname} (n={len(plot_data)}'
        if n_rejected > 0:
            title_str += f', {n_rejected} outliers)'
        else:
            title_str += ')'
        ax.set_title(title_str)
        ax.grid(True, alpha=0.3)

        # Stats overlay
        stats_text = f'μ={np.mean(plot_data):.3g}\nσ={np.std(plot_data):.3g}\nmed={np.median(plot_data):.3g}'
        ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Hide extra
    for idx in range(len(columns), len(axes_flat)):
        axes_flat[idx].axis('off')

    fig.suptitle(f'Histograms: {", ".join(columns)}', fontsize=14, fontweight='bold')
    fig.tight_layout()

    try:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f'[roman_phot] histograms -> {output_path}', file=sys.stderr)
    except Exception as exc:
        print(f'[roman_phot] WARNING: failed to save histograms: {exc}', file=sys.stderr)
    finally:
        plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description='Aperture photometry across all SCAs of a Roman WFI exposure.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        '--uri-file', metavar='PATH',
        help='Text file with one S3 (or local) URI per line (required unless --remake-mosaics)',
    )
    ap.add_argument(
        '--per-sca', action='store_true',
        help='Also write sca{NN}.csv for each SCA inside the output directory',
    )
    ap.add_argument(
        '--workers', type=int, default=8, metavar='N',
        help='Maximum parallel streaming threads (default: 8)',
    )

    # Photometry tuning (mirrors roman_view_sca)
    ap.add_argument('--fwhm',          type=float, default=1.5,  metavar='N',
                    help='Source FWHM in pixels (default: 1.5)')
    ap.add_argument('--det-sigma',     type=float, default=10.0, metavar='N',
                    help='Detection threshold in σ (default: 10.0)')
    ap.add_argument('--aper',          type=float, default=2.5,  metavar='N',
                    help='Aperture radius as multiple of FWHM (default: 2.5)')
    ap.add_argument('--annulus-inner', type=float, default=6.0,  metavar='N',
                    help='Inner annulus radius as multiple of FWHM (default: 6.0)')
    ap.add_argument('--annulus-outer', type=float, default=8.0,  metavar='N',
                    help='Outer annulus radius as multiple of FWHM (default: 8.0)')
    ap.add_argument('--bkg-poly',      type=int,   default=3,    metavar='N',
                    help='2D polynomial background degree (default: 3)')
    ap.add_argument('--snr-threshold', type=float, default=5.0,  metavar='N',
                    help='Minimum SNR to keep a source (default: 5.0)')

    # Background mosaic options
    ap.add_argument('--bkg-mosaic', action='store_true',
                    help='Produce background + source-density mosaic PNGs and save mosaic_data.npz')
    ap.add_argument('--image-mosaic', action='store_true',
                    help='Produce a focal-plane image mosaic PNG (::2 decimation, 300 dpi)')
    ap.add_argument('--remake-mosaics', metavar='PATH',
                    help='Skip photometry; regenerate mosaic PNGs from a previously saved mosaic_data.npz')
    ap.add_argument('--bkg-superpixel', type=int, default=512, metavar='N',
                    help='Superpixel bin size in pixels (default: 512)')
    ap.add_argument('--bkg-mask-sigma', type=float, default=1.5, metavar='N',
                    help='Source detection threshold for masking, in σ (default: 1.5)')
    ap.add_argument('--bkg-dilate', type=int, default=20, metavar='N',
                    help='Source mask dilation radius in pixels (default: 20)')

    # Histogram generation
    ap.add_argument('--no-hist', action='store_true',
                    help='Skip histogram generation at the end')

    args = ap.parse_args()



    # --remake-mosaics: skip photometry, load saved arrays and re-render PNGs
    if args.remake_mosaics:
        sca_maps = load_mosaic_data(args.remake_mosaics)
        out_dir = os.path.dirname(os.path.abspath(args.remake_mosaics))
        exp_label = os.path.basename(out_dir)
        make_bkg_mosaic_png(
            sca_maps, os.path.join(out_dir, 'bkg_mosaic.png'),
            superpixel=args.bkg_superpixel,
            title=f'WFI {exp_label} — background mosaic',
        )
        sources_csv = os.path.join(out_dir, 'sources.csv')
        if os.path.exists(sources_csv):
            make_source_dot_mosaic_png(
                sources_csv, os.path.join(out_dir, 'source_mosaic.png'),
                title=f'WFI {exp_label} — sources',
            )
        return

    if not args.uri_file:
        sys.exit('[roman_phot] ERROR: --uri-file is required unless --remake-mosaics is given')

    pairs = load_uri_file(args.uri_file)
    if not pairs:
        sys.exit(f'[roman_phot] ERROR: no valid URIs found in {args.uri_file}')
    print(f'[roman_phot] {len(pairs)} SCA(s) to process', file=sys.stderr)

    # Derive exposure label from the first filename: {visit_id}_{exposure_num}
    first_meta = parse_filename(pairs[0][1])
    exp_label = f'{first_meta["visit_id"]}_{first_meta["exposure_num"]}'
    exp_title = (f'WFI {first_meta["visit_id"]}  '
                 f'Exp: {first_meta["exposure_num"]}  '
                 f'{first_meta["filter"]}')
    out_dir = exp_label
    os.makedirs(out_dir, exist_ok=True)
    print(f'[roman_phot] output directory: {out_dir}/', file=sys.stderr)

    def p(name):
        """Return path inside the exposure output directory."""
        return os.path.join(out_dir, name)

    # --image-mosaic only: stream + decimate, skip all photometry
    if args.image_mosaic and not args.bkg_mosaic:
        thumbs = stream_image_mosaic(pairs, max_workers=args.workers)
        make_image_mosaic_png(
            thumbs, p('image_mosaic.png'),
            title=f'{exp_title} — image mosaic',
        )
        return

    phot_kwargs = dict(
        fwhm_pix=args.fwhm,
        detection_sigma=args.det_sigma,
        aperture_radius_fwhm=args.aper,
        annulus_inner_fwhm=args.annulus_inner,
        annulus_outer_fwhm=args.annulus_outer,
        snr_threshold=args.snr_threshold,
    )

    bkg_kwargs = None
    if args.bkg_mosaic:
        bkg_kwargs = dict(
            superpixel=args.bkg_superpixel,
            mask_sigma=args.bkg_mask_sigma,
            dilate_radius=args.bkg_dilate,
            bkg_poly_degree=args.bkg_poly,
        )

    results = phot_exposure(pairs, phot_kwargs=phot_kwargs,
                            bkg_kwargs=bkg_kwargs, out_dir=out_dir,
                            max_workers=args.workers,
                            image_mosaic=args.image_mosaic)

    if not results:
        sys.exit('[roman_phot] ERROR: no SCAs completed successfully')

    # Per-SCA CSVs
    if args.per_sca:
        for r in results:
            if r['table'] is not None and len(r['table']) > 0:
                path = p(f'sca{r["sca"]:02d}.csv')
                r['table'].write(path, format='ascii.csv', overwrite=True)
                print(f'[roman_phot] wrote {path}', file=sys.stderr)

    # Combined per-source CSV
    tables = [r['table'] for r in results if r['table'] is not None and len(r['table']) > 0]
    if tables:
        combined = vstack(tables, metadata_conflicts='silent')
        combined.write(p('sources.csv'), format='ascii.csv', overwrite=True)
        print(f'[roman_phot] combined source table ({len(combined)} rows) -> {p("sources.csv")}',
              file=sys.stderr)
    else:
        print('[roman_phot] WARNING: no sources detected in any SCA; combined CSV not written',
              file=sys.stderr)

    # Per-SCA summary CSV
    summary = build_summary(results)
    if len(summary) > 0:
        summary.write(p('summary.csv'), format='ascii.csv', overwrite=True)
        print(f'[roman_phot] summary ({len(summary)} SCAs) -> {p("summary.csv")}', file=sys.stderr)

    # Background mosaic PNG and source dot mosaic
    if args.bkg_mosaic:
        sca_maps = {r['sca']: r['bkg_map'] for r in results}
        save_mosaic_data(sca_maps, p('mosaic_data.npz'))
        make_bkg_mosaic_png(
            sca_maps, p('bkg_mosaic.png'),
            superpixel=args.bkg_superpixel,
            title=f'{exp_title} — background mosaic',
        )
        if tables:
            make_source_dot_mosaic_png(
                p('sources.csv'), p('source_mosaic.png'),
                title=f'{exp_title} — sources',
            )

    # Image mosaic
    if args.image_mosaic:
        sca_thumbs = {r['sca']: r['thumb'] for r in results}
        make_image_mosaic_png(
            sca_thumbs, p('image_mosaic.png'),
            title=f'{exp_title} — image mosaic',
        )

    # Generate histograms
    if not args.no_hist and tables:
        make_histograms_from_csv(p('sources.csv'), p('histograms.png'),
                                 columns=['aperture_sum_0', 'snr'])


if __name__ == '__main__':
    main()
