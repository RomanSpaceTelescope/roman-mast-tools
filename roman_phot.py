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

    # Custom output paths and photometry knobs:
    python roman_phot.py --uri-file my_exposure.txt \\
        --out sources.csv --summary-out summary.csv \\
        --fwhm 1.5 --det-sigma 10 --snr-threshold 5

    # Per-SCA CSVs as well:
    python roman_phot.py --uri-file my_exposure.txt --per-sca

    # Background mosaic (source-masked superpixel map in WFI focal-plane layout):
    python roman_phot.py --uri-file my_exposure.txt --bkg-mosaic
"""

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import vstack, Table

from roman_view_sca import stream_sca, run_aperture_photometry, parse_filename

# ---------------------------------------------------------------------------
# WFI focal-plane constants
# ---------------------------------------------------------------------------

# Roman WFI detector constants
_ROMAN_PIXEL_SCALE_MM = 0.01      # mm per pixel (focal-plane scale)
_ROMAN_SCA_FULL_SIZE  = 4096      # full detector size before reference-pixel removal
_ROMAN_REF_PIX        = 4         # reference pixels removed from each edge

# Focal-plane layout: SCA number → (x_center_arcmin, y_center_arcmin, rotation_deg)
# Rotation 180° means the SCA is flipped relative to the focal-plane axes.
_WFI_SCA_LAYOUT = {
     1: ( -22.14,  12.15, 180.0),
     2: ( -22.29, -37.03, 180.0),
     3: ( -22.44, -82.06,   0.0),
     4: ( -68.30,  21.84, 180.0),  # left 1.88 total, up 0.94
     5: ( -67.86, -28.28, 180.0),  # left 0.94
     6: ( -68.36, -73.06,   0.0),  # left 0.94
     7: (-112.58,  43.14, 180.0),  # left 1.88 total, up 0.94
     8: (-113.36,  -6.04, 180.0),  # left 1.88 total, up 0.94
     9: (-114.52, -50.12,   0.0),  # left 1.88 total
    10: (  22.14,  12.15, 180.0),
    11: (  22.29, -37.03, 180.0),
    12: (  22.44, -82.06,   0.0),
    13: (  66.42,  21.84, 180.0),  # right 0.94, up 0.94
    14: (  66.92, -28.28, 180.0),
    15: (  67.42, -73.06,   0.0),
    16: ( 110.70,  43.14, 180.0),  # right 0.94, up 0.94
    17: ( 111.48,  -6.98, 180.0),
    18: ( 112.64, -51.06,   0.0),
}


# ---------------------------------------------------------------------------
# Per-SCA photometry
# ---------------------------------------------------------------------------

def phot_one_sca(uri, filename, *, phot_kwargs, bkg_kwargs=None, _SourcePhotometry=None):
    """Stream one SCA and run aperture photometry.

    Returns a dict with keys: sca, detector, table, stats, sp, bkg_map.
    bkg_map is None when bkg_kwargs is None (background mosaic not requested).
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
            data, dq=dq, _SourcePhotometry=_SourcePhotometry, **phot_kwargs
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
        f'bkg={stats["bkg_level"]:.3g}  rms={stats["bkg_rms"]:.3g}  '
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

    return {
        'sca': sca_num,
        'detector': detector,
        'table': phot_table,
        'stats': stats,
        'sp': sp,
        'bkg_map': bkg_map,
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
    from astropy.stats import sigma_clipped_stats

    ny, nx = data.shape
    mask = (dq != 0) if dq is not None else np.zeros_like(data, dtype=bool)

    if data_sub is None:
        # Polynomial background fit (same logic as run_aperture_photometry)
        from astropy.modeling import models, fitting
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
        from photutils.segmentation import detect_sources
        from astropy.convolution import Gaussian2DKernel, convolve
        from scipy.ndimage import binary_dilation, generate_binary_structure, iterate_structure

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

    # Scale: focal-plane mm → canvas superpixels
    # 1 superpixel = superpixel * _ROMAN_PIXEL_SCALE_MM mm
    sp_mm = superpixel * _ROMAN_PIXEL_SCALE_MM   # mm per superpixel (5.12 mm for 512 px)
    scale = 1.0 / sp_mm                          # superpixels per mm

    # Infer tile shape from available maps
    tile_shapes = [v.shape for v in sca_maps.values() if v is not None]
    if not tile_shapes:
        print(f'[roman_phot] WARNING: no background maps to render; skipping PNG', file=sys.stderr)
        return
    tile_h, tile_w = tile_shapes[0]

    # Canvas bounds in arcminutes (with half-tile padding)
    coords = list(_WFI_SCA_LAYOUT.values())
    all_cx, all_cy = zip(*[(cx, cy) for cx, cy, _ in coords])
    half_w = (tile_w / scale) / 2
    half_h = (tile_h / scale) / 2
    x_min = min(all_cx) - half_w
    x_max = max(all_cx) + half_w
    y_min = min(all_cy) - half_h
    y_max = max(all_cy) + half_h

    canvas_w = int(np.ceil((x_max - x_min) * scale)) + 2
    canvas_h = int(np.ceil((y_max - y_min) * scale)) + 2
    canvas = np.full((canvas_h, canvas_w), np.nan)

    for sca_num, (cx_am, cy_am, rot) in _WFI_SCA_LAYOUT.items():
        tile = sca_maps.get(sca_num)
        if tile is None:
            continue

        if rot == 180.0:
            tile = np.rot90(tile, 2)

        # Convert centre from arcmin to canvas pixel (y flipped: sky up → image down)
        cx_px = (cx_am - x_min) * scale
        cy_px = (y_max - cy_am) * scale

        # Position tile so its center is at (cx_px, cy_px). Keep float precision
        # until placement to avoid accumulating rounding errors.
        col0_float = cx_px - tile_w / 2
        row0_float = cy_px - tile_h / 2
        col0 = int(np.floor(col0_float))
        row0 = int(np.floor(row0_float))
        col1 = col0 + tile_w
        row1 = row0 + tile_h

        # Clip to canvas
        r0 = max(row0, 0); r1 = min(row1, canvas_h)
        c0 = max(col0, 0); c1 = min(col1, canvas_w)
        tr0 = r0 - row0; tc0 = c0 - col0
        canvas[r0:r1, c0:c1] = tile[tr0:tr0 + (r1 - r0), tc0:tc0 + (c1 - c0)]

    # Colour stretch: symmetric around zero (diverging)
    finite = canvas[np.isfinite(canvas)]
    if len(finite) == 0:
        print(f'[roman_phot] WARNING: background mosaic is entirely NaN; skipping PNG', file=sys.stderr)
        return
    abs_lim = max(abs(np.percentile(finite, pct_lo)), abs(np.percentile(finite, pct_hi)))
    vmin, vmax = -abs_lim, abs_lim

    fig, ax = plt.subplots(figsize=(14, 8), facecolor='#1a1a1a')
    ax.set_facecolor('#1a1a1a')

    im = ax.imshow(
        canvas, origin='upper', cmap='RdBu_r',
        vmin=vmin, vmax=vmax, interpolation='nearest',
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label('Background residual (DN/s)', color='white', fontsize=10)
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')

    # Label each SCA
    for sca_num, (cx_am, cy_am, _) in _WFI_SCA_LAYOUT.items():
        cx_px = (cx_am - x_min) * scale
        cy_px = (y_max - cy_am) * scale
        ax.text(cx_px, cy_px, f'{sca_num:02d}',
                ha='center', va='center', fontsize=7,
                color='white', alpha=0.6, fontweight='bold')

    ax.set_title(title or 'Roman WFI — source-masked background mosaic',
                 color='white', fontsize=12, pad=10)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'[roman_phot] background mosaic -> {out_path}', file=sys.stderr)


# ---------------------------------------------------------------------------
# Parallel fan-out
# ---------------------------------------------------------------------------

def phot_exposure(uri_filename_pairs, *, phot_kwargs=None, bkg_kwargs=None, max_workers=8):
    """Run photometry on a list of (uri, filename) pairs in parallel.

    Returns a list of result dicts (sorted by SCA number), with None entries
    removed.
    """
    # Load roman_lolo once at the start, shared across all SCAs
    try:
        from roman_lolo.romanphot import SourcePhotometry
    except ImportError:
        raise ImportError('roman-lolo not installed; install with: pip install -e .')

    if phot_kwargs is None:
        phot_kwargs = {}

    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(phot_one_sca, uri, fn,
                        phot_kwargs=phot_kwargs, bkg_kwargs=bkg_kwargs,
                        _SourcePhotometry=SourcePhotometry): fn
            for uri, fn in uri_filename_pairs
        }
        for fut in as_completed(futures):
            res = fut.result()
            if res is not None:
                results[res['sca']] = res

    return [results[k] for k in sorted(results)]


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
    from astropy.table import Table

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
            'bkg_rms':           stats['bkg_rms'],
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
        '--uri-file', required=True, metavar='PATH',
        help='Text file with one S3 (or local) URI per line',
    )
    ap.add_argument(
        '--out', default='roman_phot_combined.csv', metavar='PATH',
        help='Combined per-source CSV output (default: roman_phot_combined.csv)',
    )
    ap.add_argument(
        '--summary-out', default='roman_phot_summary.csv', metavar='PATH',
        help='Per-SCA summary CSV output (default: roman_phot_summary.csv)',
    )
    ap.add_argument(
        '--per-sca', action='store_true',
        help='Also write roman_phot_sca{NN}.csv for each SCA',
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
                    help='Produce a source-masked superpixel background mosaic PNG')
    ap.add_argument('--bkg-mosaic-out', default='roman_phot_bkg_mosaic.png', metavar='PATH',
                    help='Background mosaic PNG output (default: roman_phot_bkg_mosaic.png)')
    ap.add_argument('--bkg-superpixel', type=int, default=512, metavar='N',
                    help='Superpixel bin size in pixels (default: 512)')
    ap.add_argument('--bkg-mask-sigma', type=float, default=1.5, metavar='N',
                    help='Source detection threshold for masking, in σ (default: 1.5)')
    ap.add_argument('--bkg-dilate', type=int, default=20, metavar='N',
                    help='Source mask dilation radius in pixels (default: 20)')

    # Histogram generation
    ap.add_argument('--no-hist', action='store_true',
                    help='Skip histogram generation at the end')
    ap.add_argument('--hist-output', metavar='PATH',
                    help='Histogram output path (default: roman_phot_histograms.png)')

    args = ap.parse_args()

    pairs = load_uri_file(args.uri_file)
    if not pairs:
        sys.exit(f'[roman_phot] ERROR: no valid URIs found in {args.uri_file}')
    print(f'[roman_phot] {len(pairs)} SCA(s) to process', file=sys.stderr)

    phot_kwargs = dict(
        fwhm_pix=args.fwhm,
        detection_sigma=args.det_sigma,
        aperture_radius_fwhm=args.aper,
        annulus_inner_fwhm=args.annulus_inner,
        annulus_outer_fwhm=args.annulus_outer,
        bkg_poly_degree=args.bkg_poly,
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
                            bkg_kwargs=bkg_kwargs, max_workers=args.workers)

    if not results:
        sys.exit('[roman_phot] ERROR: no SCAs completed successfully')

    # Per-SCA CSVs
    if args.per_sca:
        for r in results:
            if r['table'] is not None and len(r['table']) > 0:
                path = f'roman_phot_sca{r["sca"]:02d}.csv'
                r['table'].write(path, format='ascii.csv', overwrite=True)
                print(f'[roman_phot] wrote {path}', file=sys.stderr)

    # Combined per-source CSV
    tables = [r['table'] for r in results if r['table'] is not None and len(r['table']) > 0]
    if tables:
        combined = vstack(tables)
        combined.write(args.out, format='ascii.csv', overwrite=True)
        print(f'[roman_phot] combined source table ({len(combined)} rows) -> {args.out}',
              file=sys.stderr)
    else:
        print('[roman_phot] WARNING: no sources detected in any SCA; combined CSV not written',
              file=sys.stderr)

    # Per-SCA summary CSV
    summary = build_summary(results)
    if len(summary) > 0:
        summary.write(args.summary_out, format='ascii.csv', overwrite=True)
        print(f'[roman_phot] summary ({len(summary)} SCAs) -> {args.summary_out}', file=sys.stderr)

    # Background mosaic PNG
    if args.bkg_mosaic:
        sca_maps = {r['sca']: r['bkg_map'] for r in results}
        import os
        exp_label = os.path.splitext(os.path.basename(args.uri_file))[0]
        make_bkg_mosaic_png(
            sca_maps, args.bkg_mosaic_out,
            superpixel=args.bkg_superpixel,
            title=f'Roman WFI background mosaic — {exp_label}',
        )

    # Generate histograms
    if not args.no_hist and tables:
        hist_output = args.hist_output or 'roman_phot_histograms.png'
        make_histograms_from_csv(args.out, hist_output, columns=['aperture_sum_0', 'snr'])


if __name__ == '__main__':
    main()
