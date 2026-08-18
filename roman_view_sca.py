#!/usr/bin/env python
"""Stream a single Roman WFI SCA from MAST or S3 into DS9 or matplotlib with channel grid overlay.

Optionally runs aperture photometry and overlays detected sources with SNR-coloured apertures.

DS9 must already be running before you call this script (for --display ds9).

Usage
-----
    # Via roman_mast query (requires MAST auth token):
    # List matching exposures:
    conda run -n roman-mast-tools python roman_view_sca.py --program 114 --pass 57 --list

    # View SCA 1 from first matching exposure in DS9:
    conda run -n roman-mast-tools python roman_view_sca.py --program 114 --pass 57 --sca 1

    # View SCA 11 from 2nd matching exposure with photometry:
    conda run -n roman-mast-tools python roman_view_sca.py \\
        --program 114 --pass 57 --exposure 2 --sca 11 --phot --phot-out phot.csv

    # Direct S3 streaming (anonymous access):
    python roman_view_sca.py \\
        s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \\
        r0003201001001001004_0001_wfi11_f106_cal.asdf

    # Anonymous S3 with matplotlib + photometry:
    python roman_view_sca.py --display mpl --phot --phot-out phot.csv \\
        s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \\
        r0003201001001001004_0001_wfi11_f106_cal.asdf

    # Suppress channel dividers, keep section grid:
    python roman_view_sca.py s3://... r0003...asdf --no-channels

    # Connect to a named DS9 instance:
    python roman_view_sca.py --program 114 --pass 57 --sca 1 --ds9 myds9

    # Background subtraction only (no photometry) — shows model + residuals figure:
    python roman_view_sca.py --display mpl --bkg \\
        s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \\
        r0003201001001001004_0001_wfi11_f106_cal.asdf

    # Background + photometry together (shares the same polynomial fit):
    python roman_view_sca.py --display mpl --bkg --phot \\
        s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \\
        r0003201001001001004_0001_wfi11_f106_cal.asdf
"""

import argparse
import io
import os
import sys

import numpy as np
import s3fs
from astropy.io import fits
from astropy.wcs import WCS

# ---------------------------------------------------------------------------
# Focal-plane rotation table (from MPA_SCA_info).
# SCAs 3, 6, 9, 12, 15, 18 are not rotated (r=0); all others are r=180.
# ---------------------------------------------------------------------------
_SCA_ROTATION = {n: (0 if n % 3 == 0 else 180) for n in range(1, 19)}


# ---------------------------------------------------------------------------
# Streaming
# ---------------------------------------------------------------------------

def stream_sca(uri, filename=None, *, sip_degree=4, asdf_file=None):
    """Open an ASDF file and return (data_f32, detector, wcs_hdr, dq).

    Three input modes:
    1. Via open AsdfFile (asdf_file): datamodel already loaded, no network I/O.
    2. Via S3 URI + filename: 's3://...' → anonymous access.
    3. Via local path + filename: '/path/to/dir' + 'filename.asdf'.

    Materialises ``dm.data`` and computes the SIP WCS approximation while the
    ASDF file handle is still open (gwcs needs the live tree for to_fits_sip).
    """
    import roman_datamodels as rdm

    if asdf_file is not None:
        # Already have an open AsdfFile (from roman_mast streaming)
        print(f'[view_sca] using provided AsdfFile', file=sys.stderr)
        dm = rdm.open(asdf_file)
    elif uri.startswith('s3://'):
        # S3 file path (anonymous access)
        path = uri.rstrip('/') + '/' + filename
        fs = s3fs.S3FileSystem(anon=True)
        print(f'[view_sca] streaming from S3: {path}', file=sys.stderr)
        with fs.open(path, 'rb') as f:
            dm = rdm.open(f)
            data     = np.asarray(dm.data[...], dtype=np.float32)
            detector = str(dm.meta.instrument.detector).upper()
            try:
                dq = np.asarray(dm.dq[...], dtype=np.int32)
            except AttributeError:
                dq = None   # uncal products have no DQ layer
            print(f'[view_sca] computing SIP WCS (degree={sip_degree}) ...', file=sys.stderr)
            gwcs    = dm.meta.wcs
            wcs_hdr = gwcs.to_fits_sip(bounding_box=gwcs.bounding_box, degree=sip_degree)
        print(f'[view_sca] {detector}  shape={data.shape}', file=sys.stderr)
        return data, detector, wcs_hdr, dq
    else:
        # Local file path — materialise inside the with block; ASDF blocks are
        # lazily mapped and become inaccessible once the file handle is closed.
        path = os.path.join(uri.rstrip('/'), filename)
        print(f'[view_sca] reading from local: {path}', file=sys.stderr)
        with open(path, 'rb') as f:
            dm = rdm.open(f)
            data     = np.asarray(dm.data[...], dtype=np.float32)
            detector = str(dm.meta.instrument.detector).upper()
            try:
                dq = np.asarray(dm.dq[...], dtype=np.int32)
            except AttributeError:
                dq = None
            print(f'[view_sca] computing SIP WCS (degree={sip_degree}) ...', file=sys.stderr)
            gwcs    = dm.meta.wcs
            wcs_hdr = gwcs.to_fits_sip(bounding_box=gwcs.bounding_box, degree=sip_degree)
        print(f'[view_sca] {detector}  shape={data.shape}', file=sys.stderr)
        return data, detector, wcs_hdr, dq

    # asdf_file path (already-open handle from roman_mast streaming)
    data     = np.asarray(dm.data[...], dtype=np.float32)
    detector = str(dm.meta.instrument.detector).upper()
    try:
        dq = np.asarray(dm.dq[...], dtype=np.int32)
    except AttributeError:
        dq = None
    print(f'[view_sca] computing SIP WCS (degree={sip_degree}) ...', file=sys.stderr)
    gwcs    = dm.meta.wcs
    wcs_hdr = gwcs.to_fits_sip(bounding_box=gwcs.bounding_box, degree=sip_degree)
    print(f'[view_sca] {detector}  shape={data.shape}', file=sys.stderr)
    return data, detector, wcs_hdr, dq


# ---------------------------------------------------------------------------
# Filename parsing
# ---------------------------------------------------------------------------

def parse_filename(filename):
    """Extract visit ID, exposure number, SCA number, and filter from Roman filename.

    E.g. r0003201001001001004_0001_wfi11_f106_cal.asdf
    """
    base = filename.replace('.asdf', '')
    parts = base.split('_')

    visit_id = parts[0]           # r0003201001001001004
    exposure_num = parts[1]       # 0001
    sca_str = parts[2]            # wfi11
    filter_str = parts[3]         # f106

    sca_num = int(sca_str.replace('wfi', ''))

    return {
        'visit_id': visit_id,
        'exposure_num': exposure_num,
        'sca_num': sca_num,
        'filter': filter_str,
    }


# ---------------------------------------------------------------------------
# Region string builder
# ---------------------------------------------------------------------------

def make_channel_regions(data, detector, *, show_channels=False, show_gridlines=False):
    """Return a DS9 region string (image coordinates, 1-indexed) for the
    H4RG readout channel and section grid.

    Parameters
    ----------
    data : ndarray (ny, nx)
        The already-materialised science array.  Used only for shape.
    detector : str
        e.g. ``'WFI11'``.
    show_channels : bool
        Draw 32 thin channel dividers (128 px each) and channel numbers.
    show_gridlines : bool
        Draw 8-section major grid (512 px each) and section labels a–h / A–H.
    """
    ny, nx   = data.shape
    ref_pix  = (4096 - nx) // 2   # 4 if reference pixels stripped, else 0
    sca_num  = int(detector[-2:])
    rotated  = _SCA_ROTATION[sca_num] == 180

    rgns = ['# Region file format: DS9 version 4.1', 'image']

    # 32 readout channel dividers — thin, dashed
    if show_channels:
        for j in range(1, 32):
            x = 128 * j - ref_pix + 1          # 1-indexed pixel boundary
            rgns.append(
                f'line({x},1,{x},{ny})'
                f' # color=white width=1 line=0 0 dash=1'
            )

    section_xs = [512 * j - ref_pix + 1 for j in range(1, 8)]

    # Major section grid — thicker, solid
    if show_gridlines:
        for x in section_xs:
            rgns.append(f'line({x},1,{x},{ny}) # color=white width=2 line=0 0')
        for j in range(1, 8):
            y = 512 * j - ref_pix + 1
            rgns.append(f'line(1,{y},{nx},{y}) # color=white width=2 line=0 0')

        # Section labels: a–h right of image, A–H above
        for i, (rl, cl) in enumerate(zip('abcdefgh', 'ABCDEFGH')):
            mid = 256 + i * 512 - ref_pix + 1
            rgns.append(
                f'text({nx + 60},{mid})'
                f' # text={{{rl}}} color=black font="helvetica 10 bold roman"'
            )
            rgns.append(
                f'text({mid},{ny + 60})'
                f' # text={{{cl}}} color=black font="helvetica 10 bold roman"'
            )

    # Channel numbers 1–32 (or 32–1 for r=0 SCAs) near the top edge
    if show_channels:
        for i in range(32):
            chan_num = i + 1 if rotated else 32 - i
            x = 64 + i * 128 - ref_pix + 1
            rgns.append(
                f'text({x},{ny + 110})'
                f' # text={{{chan_num}}} color=black font="helvetica 7 normal roman"'
            )

    return '\n'.join(rgns)




def phot_to_region_str(sources, phot_table, sp):
    """Generate a DS9 region string for aperture photometry results.

    DAOStarFinder returns 0-based pixel coords (numpy convention, center of first
    pixel = (0, 0)). DS9 'image' coordinates are FITS 1-based (center of first
    pixel = (1, 1)), so we add 1.0 for correct alignment.

    Returns a multi-line region string with aperture circles SNR-colour-coded
    (no labels; see CSV for details).
    """
    rgns = ['# Region file format: DS9 version 4.1', 'image']

    for i in range(len(sources)):
        x = sources['x_centroid'][i] + 1.0   # numpy 0-based -> DS9 image (FITS) 1-based
        y = sources['y_centroid'][i] + 1.0

        # SNR-based colour
        color = 'red'
        if phot_table is not None and i < len(phot_table) and 'snr' in phot_table.colnames:
            snr_val = phot_table['snr'][i]
            if snr_val >= 50:
                color = 'green'
            elif snr_val >= 20:
                color = 'cyan'
            elif snr_val >= 10:
                color = 'yellow'
            elif snr_val >= 5:
                color = 'magenta'

        # Aperture circle (no text label)
        rgns.append(
            f'circle({x:.3f},{y:.3f},{sp.aperture_radius:.2f}) '
            f'# color={color} width=2'
        )

        # Background annulus (dashed)
        rgns.append(
            f'annulus({x:.3f},{y:.3f},'
            f'{sp.annulus_inner:.2f},{sp.annulus_outer:.2f}) '
            f'# color={color} width=1 dash=1'
        )

    return '\n'.join(rgns)


# ---------------------------------------------------------------------------
# DS9 display
# ---------------------------------------------------------------------------

def display_in_ds9(data, regions, *, title='Roman WFI', wcs_header=None,
                   dq=None, ds9_target=None, scale='zscale', cmap='viridis',
                   sources=None, phot_table=None, sp=None, data_sub=None):
    """Send a 2-D array to DS9 and apply the region overlay.

    Parameters
    ----------
    data : ndarray (ny, nx), float32
    regions : str
        DS9 region string from ``make_channel_regions``.
    title : str
        Value written to the FITS OBJECT keyword (shown in DS9 title bar).
    ds9_target : str, optional
        XPA target name of a running DS9 instance.  None → default instance.
    scale : str
        DS9 scale algorithm (``'zscale'``, ``'log'``, ``'linear'``, …).
    cmap : str
        DS9 colour map name.
    """
    try:
        import pyds9
    except ImportError:
        sys.exit(
            'pyds9 is not installed in this environment.\n'
            'Install with: pip install pyds9'
        )

    try:
        d = pyds9.DS9(target=ds9_target) if ds9_target else pyds9.DS9()
    except Exception as e:
        sys.exit(
            f'Could not connect to DS9: {e}\n'
            'Make sure DS9 is running first:  ds9 &'
        )

    # Build header: include SIP WCS for RA/Dec display, but also keep NAXIS/CRPIX for region alignment
    hdr = wcs_header.copy() if wcs_header is not None else fits.Header()
    hdr['OBJECT'] = title
    buf = io.BytesIO()
    fits.PrimaryHDU(data=data, header=hdr).writeto(buf)
    fits_bytes = buf.getvalue()

    print(f'[view_sca_ds9] sending {len(fits_bytes)/1e6:.1f} MB to DS9', file=sys.stderr)
    d.set('frame delete all')
    d.set('frame new')
    d.set('fits', fits_bytes)
    d.set(f'scale {scale}')
    d.set(f'cmap {cmap}')
    d.set('zoom to fit')

    if dq is not None:
        dq_buf = io.BytesIO()
        fits.PrimaryHDU(data=dq).writeto(dq_buf)
        d.set('mask clear')
        d.set('mask color red')
        d.set('mask transparency 50')
        d.set('mask mark nonzero')
        d.set('fits mask', dq_buf.getvalue())

    d.set('regions delete all')
    d.set('regions', regions)

    if sources is not None and len(sources) > 0 and sp is not None:
        phot_regions = phot_to_region_str(sources, phot_table, sp)
        # Append photometry regions to existing channel/grid regions
        combined_regions = regions + '\n' + phot_regions
        d.set('regions delete all')
        d.set('regions', combined_regions)

    if data_sub is not None:
        buf2 = io.BytesIO()
        hdr2 = wcs_header.copy() if wcs_header is not None else fits.Header()
        hdr2['OBJECT'] = title + ' [bkg residuals]'
        fits.PrimaryHDU(data=data_sub, header=hdr2).writeto(buf2)
        d.set('frame new')
        d.set('fits', buf2.getvalue())
        d.set('scale zscale')
        d.set('cmap bb')
        d.set('zoom to fit')
        print('[view_sca_ds9] residuals pushed to frame 2 (blink with frame 1)', file=sys.stderr)
        d.set('frame first')

    print('[view_sca_ds9] done — DS9 is open', file=sys.stderr)
    return d


# ---------------------------------------------------------------------------
# Matplotlib display
# ---------------------------------------------------------------------------

def display_in_mpl(data, detector, *, dq=None, title=None, wcs_header=None,
                   show_channels=False, show_gridlines=False,
                   figsize=(9, 8), norm=None, sources=None, phot_table=None, sp=None, stats=None,
                   save_path=None):
    """Display a single SCA in an interactive matplotlib window.

    When wcs_header is provided the axes use WCSAxes projection so ticks and
    the status bar show RA/Dec; the data value (and DQ flag if present) are
    appended to the WCSAxes format_coord output.
    """
    from astropy.visualization import simple_norm

    ny, nx  = data.shape
    ref_pix = (4096 - nx) // 2
    sca_num = int(detector[-2:])
    rotated = _SCA_ROTATION[sca_num] == 180

    if norm is None:
        norm = simple_norm(data, 'asinh', vmin=0.5, vmax=4)

    projection = WCS(wcs_header) if wcs_header is not None else None

    if save_path is not None:
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        fig = Figure(figsize=figsize)
        FigureCanvasAgg(fig)
        ax = fig.add_subplot(111, projection=projection)
    else:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=figsize,
                               subplot_kw={'projection': projection} if projection else {})
    ax.imshow(data, norm=norm, origin='lower')
    if dq is not None:
        dq_overlay = np.where(dq != 0, 1.0, np.nan)
        ax.imshow(dq_overlay, origin='lower', cmap='Reds', alpha=0.5, vmin=0, vmax=1)

    ax.set(
        xlabel='X Science Axis (pixels)',
        ylabel='Y Science Axis (pixels)',
        title=title or f'Roman WFI {detector}',
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    if hasattr(ax, 'coords'):          # WCSAxes draws its own frame on top of spines
        ax.coords.frame.set_linewidth(0)

    if show_channels:
        for j in range(1, 32):
            ax.axvline(128 * j - ref_pix, color='white', alpha=0.2, lw=0.5)

    section_locs = [512 * j - ref_pix for j in range(1, 8)]
    if show_gridlines:
        for pos in section_locs:
            ax.axvline(pos, color='white', alpha=0.4, lw=0.8)
            ax.axhline(pos, color='white', alpha=0.4, lw=0.8)

    if projection is None:
        # Plain axes: show native pixel coordinates at section boundaries
        # Disable tick labels when saving to PNG to avoid expensive text measurement in Agg backend
        if save_path is None:
            tick_locs = [0] + section_locs + [nx - 1]
            tick_vals = [ref_pix] + [512 * j for j in range(1, 8)] + [4096 - ref_pix]
            ax.set_xticks(tick_locs)
            ax.set_yticks(tick_locs)
            if rotated:
                ax.set_xticklabels(tick_vals)
                ax.set_yticklabels(tick_vals[::-1])
            else:
                ax.set_xticklabels(tick_vals[::-1])
                ax.set_yticklabels(tick_vals)
        else:
            ax.set_xticks([])
            ax.set_yticks([])

    col_label_y = ny + 20
    ax.set_xlim(-0.5, nx + 50)
    ax.set_ylim(-0.5, ny + (120 if show_channels else 60))

    if show_gridlines:
        for i, (rl, cl) in enumerate(zip('abcdefgh', 'ABCDEFGH')):
            mid = 256 + i * 512 - ref_pix
            ax.text(nx + 18, mid, rl, va='center', fontsize=9)
            ax.text(mid, col_label_y, cl, ha='center', fontsize=9)

    if show_channels:
        for i in range(32):
            chan_num = i + 1 if rotated else 32 - i
            ax.text(64 + i * 128 - ref_pix, ny + 100, str(chan_num),
                    ha='center', fontsize=6)

    if sources is not None and len(sources) > 0 and sp is not None:
        from matplotlib.collections import EllipseCollection

        xs = np.asarray(sources['x_centroid'], dtype=float)
        ys = np.asarray(sources['y_centroid'], dtype=float)

        # Assign SNR-based colours per source
        colors = np.full(len(sources), 'red', dtype=object)
        if phot_table is not None and 'snr' in phot_table.colnames:
            snr = np.asarray(phot_table['snr'][:len(sources)], dtype=float)
            colors[snr >= 5]  = 'magenta'
            colors[snr >= 10] = 'yellow'
            colors[snr >= 20] = 'cyan'
            colors[snr >= 50] = 'green'

        offsets = np.column_stack([xs, ys])
        ax.set_autoscale_on(False)

        for color in np.unique(colors):
            mask = colors == color
            off = offsets[mask]
            n = mask.sum()
            for radius, lw, ls, alpha in [
                (sp.aperture_radius, 1.5, '-',  0.8),
                (sp.annulus_inner,   1.0, '--', 0.4),
                (sp.annulus_outer,   1.0, '--', 0.4),
            ]:
                col = EllipseCollection(
                    widths=2*radius, heights=2*radius, angles=0, units='xy',
                    offsets=off, offset_transform=ax.transData,
                    facecolors='none', edgecolors=color,
                    linewidths=lw, linestyles=ls, alpha=alpha,
                )
                ax.add_collection(col)

        ax.set_autoscale_on(True)

    if stats is not None and save_path is None:
        # Display background statistics in a text panel (upper-left corner)
        # Skip for saved PNGs to avoid expensive text measurement
        stats_text = (
            f"Background Stats:\n"
            f"Level: {stats['bkg_level']:.2f}\n"
            f"RMS: {stats.get('bkg_rms_median', np.median(stats['bkg_rms'])):.2f}\n"
            f"Threshold: {stats['threshold']:.2f}\n"
            f"Sources: {stats['n_sources']}"
        )
        ax.text(0.02, 0.98, stats_text,
               transform=ax.transAxes,
               verticalalignment='top', horizontalalignment='left',
               fontsize=9, family='monospace',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    _wcs_fmt = ax.format_coord   # WCSAxes provides RA/Dec; plain axes gives x/y
    def _fmt(x, y):
        col, row = int(round(x)), int(round(y))
        # WCSAxes can raise when the WCS transform fails near image edges
        try:
            wcs_str = _wcs_fmt(x, y)
        except Exception:
            wcs_str = f'x={x:.1f}  y={y:.1f}'
        if 0 <= row < ny and 0 <= col < nx:
            val = data[row, col]
            dq_str = f'  dq={dq[row, col]}' if dq is not None else ''
            return wcs_str + f'  x={col + ref_pix}  y={row + ref_pix}  val={val:.4g}{dq_str}'
        return wcs_str
    ax.format_coord = _fmt

    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=300)
    else:
        plt.show()
    return fig, ax


def display_residuals_mpl(data_sub, detector, *, bkg_fit=None, bkg_level=0.0,
                           residual_rms=None, poly_degree=None, dq=None,
                           wcs_header=None, figsize=(14, 6), n_sigma=5.0,
                           save_path=None):
    """Two-panel diagnostic figure for background subtraction quality.

    Left panel  — polynomial background model (viridis, full range).
    Right panel — data minus background (RdBu_r, symmetric ±n_sigma × RMS).

    Pass bkg_fit=None to show only the residuals panel.
    """
    import matplotlib.colors as mcolors

    if residual_rms is None:
        sky_mask = np.isfinite(data_sub)
        if dq is not None:
            sky_mask &= (dq == 0)
        flat = data_sub[sky_mask]
        if flat.size == 0:
            flat = data_sub[np.isfinite(data_sub)].ravel()
        _, _, residual_rms = sigma_clipped_stats(flat, sigma=5.0)

    half = max(n_sigma * residual_rms, 1e-10)
    resid_norm = mcolors.Normalize(vmin=-half, vmax=half)

    n_panels = 2 if bkg_fit is not None else 1

    if save_path is not None:
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        fig = Figure(figsize=figsize)
        FigureCanvasAgg(fig)
        axes = [fig.add_subplot(1, n_panels, i + 1) for i in range(n_panels)]
    else:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, n_panels, figsize=figsize)
        if n_panels == 1:
            axes = [axes]

    if bkg_fit is not None:
        im0 = axes[0].imshow(bkg_fit, origin='lower', cmap='viridis',
                             aspect='equal', interpolation='nearest')
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04,
                     label='Background (DN/s)')
        deg_str = f'  degree={poly_degree}' if poly_degree is not None else ''
        axes[0].set_title(f'Background model{deg_str}')
        axes[0].set_xlabel('X (pixels)')
        axes[0].set_ylabel('Y (pixels)')

    im1 = axes[-1].imshow(data_sub, origin='lower', cmap='RdBu_r',
                          norm=resid_norm, aspect='equal', interpolation='nearest')
    fig.colorbar(im1, ax=axes[-1], fraction=0.046, pad=0.04,
                 label='Residual (DN/s)')
    axes[-1].set_title(
        f'Residuals after background subtraction\n'
        f'RMS = {residual_rms:.4g} DN/s    '
        f'level = {bkg_level:.4g} DN/s    '
        f'stretch = ±{n_sigma:.0f}σ = ±{half:.3g}'
    )
    axes[-1].set_xlabel('X (pixels)')
    if n_panels == 1:
        axes[-1].set_ylabel('Y (pixels)')

    if dq is not None:
        dq_overlay = np.where(dq != 0, 1.0, np.nan)
        for ax in axes:
            ax.imshow(dq_overlay, origin='lower', cmap='Reds',
                      alpha=0.35, vmin=0, vmax=1, aspect='equal')

    fig.suptitle(f'Roman WFI {detector} — background subtraction quality',
                 fontsize=12)
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=300)
    else:
        plt.show()
    return fig, axes


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=(
            'Stream a single Roman WFI SCA from MAST (via roman_mast) or S3 and display it '
            'in DS9 or matplotlib with H4RG channel/section grid overlay.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Roman MAST query mode (optional positional args for S3/local fallback)
    ap.add_argument('uri', nargs='?', default=None,
                    help='S3 base URI (e.g. s3://stpubdata/roman/...) or local directory path (e.g. ./cache)')
    ap.add_argument('filename', nargs='?', default=None,
                    help='ASDF filename, e.g. r0003..._wfi11_f106_cal.asdf')

    # Roman MAST filters (when not using uri + filename)
    ap.add_argument('--program', type=int, default=None, metavar='N',
                    help='APT program ID (e.g. 114)')
    ap.add_argument('--pass', type=int, default=None, metavar='N', dest='pass_',
                    help='Pass within the program (e.g. 57)')
    ap.add_argument('--execution-plan', type=int, default=None, metavar='N',
                    help='Execution plan within the program')
    ap.add_argument('--segment', type=int, default=None, metavar='N',
                    help='Segment within the pass')
    ap.add_argument('--observation', type=int, default=None, metavar='N',
                    help='Observation within the segment')
    ap.add_argument('--visit', type=int, default=None, metavar='N',
                    help='Visit within the observation')
    ap.add_argument('--exposure', type=int, default=1, metavar='N',
                    help='Exposure index (1-based, matching --list order; default: 1)')
    ap.add_argument('--sca', type=int, required=False, metavar='N',
                    help='SCA number to view (1-18). Required when using roman_mast query mode.')
    ap.add_argument('--optical-element', default=None, metavar='ELEM',
                    help='Filter name (e.g. F106, F129)')
    ap.add_argument('--detector', default=None, metavar='NAME',
                    help='Detector name (e.g. WFI04)')
    ap.add_argument('--data-level', type=int, choices=[1, 2], default=2, metavar='N',
                    help='Data level: 1=uncal, 2=cal (default: 2)')
    ap.add_argument('--list', action='store_true', dest='list_exposures',
                    help='List matching exposures and exit (no display)')

    # Display options (all modes)
    ap.add_argument('--channels', action='store_true',
                    help='Draw 32-channel thin dividers and channel numbers')
    ap.add_argument('--grid',     action='store_true',
                    help='Draw 8-section major grid and section labels')
    ap.add_argument('--display', choices=['ds9', 'mpl'], default='ds9',
                    help='Display backend (default: ds9)')
    ap.add_argument('--no-dq',        action='store_true',
                    help='Skip the bad-pixel mask overlay')
    ap.add_argument('--sip-degree', type=int, default=4, metavar='N',
                    help='SIP polynomial degree for gwcs->FITS WCS (default: 4; ds9 only)')
    ap.add_argument('--scale',   default='zscale',
                    help='DS9 scale algorithm (default: zscale; ds9 only)')
    ap.add_argument('--cmap',    default='viridis',
                    help='DS9 colour map (default: viridis; ds9 only)')
    ap.add_argument('--ds9',     default=None, metavar='TARGET',
                    help='XPA target name of a running DS9 (default: any; ds9 only)')

    # Background subtraction options
    ap.add_argument('--bkg', action='store_true',
                    help='Fit and subtract a 2D polynomial background; show a two-panel '
                         'residuals figure (model + residuals). Also usable with --phot '
                         '(shares the same polynomial fit).')
    ap.add_argument('--bkg-scale', type=float, default=5.0, metavar='N',
                    help='Residuals plot colour stretch in units of sigma-clipped RMS '
                         '(default: 5.0; --bkg only)')

    # Photometry options
    ap.add_argument('--phot', action='store_true',
                    help='Enable aperture photometry analysis')
    ap.add_argument('--fwhm', type=float, default=1.5, metavar='N',
                    help='Source FWHM in pixels (default: 1.5; photometry only)')
    ap.add_argument('--det-sigma', type=float, default=10.0, metavar='N',
                    help='DAOStarFinder detection threshold σ (default: 10.0; photometry only)')
    ap.add_argument('--aper', type=float, default=2.5, metavar='N',
                    help='Aperture radius as multiple of FWHM (default: 2.5; photometry only)')
    ap.add_argument('--annulus-inner', type=float, default=6.0, metavar='N',
                    help='Inner annulus radius as multiple of FWHM (default: 6.0; photometry only)')
    ap.add_argument('--annulus-outer', type=float, default=8.0, metavar='N',
                    help='Outer annulus radius as multiple of FWHM (default: 8.0; photometry only)')
    ap.add_argument('--bkg-poly-degree', type=int, default=3, metavar='N',
                    help='2D polynomial background fitting degree (default: 3; --bkg and --phot)')
    ap.add_argument('--snr-threshold', type=float, default=5.0, metavar='N',
                    help='Minimum SNR to include a source in output (default: 5.0; photometry only)')
    ap.add_argument('--phot-out', default=None, metavar='PATH',
                    help='Write photometry table to this CSV path (photometry only)')

    args = ap.parse_args()

    # Determine which mode: roman_mast query or direct S3/local
    if args.program is not None or args.pass_ is not None or args.execution_plan is not None or \
       args.segment is not None or args.observation is not None or args.visit is not None or \
       args.detector is not None or args.optical_element is not None:
        # Roman MAST mode
        from roman_mast import list_data, print_summary

        print(f'[view_sca] Querying MAST...', file=sys.stderr)
        res = list_data(
            program=args.program,
            execution_plan=args.execution_plan,
            pass_=args.pass_,
            segment=args.segment,
            observation=args.observation,
            visit=args.visit,
            detector=args.detector,
            optical_element=args.optical_element,
            sca_only=True,
            data_level=args.data_level,
        )
        if len(res.exposures) == 0:
            sys.exit('[view_sca] No exposures found matching your query')

        # List mode: show matching exposures and exit
        if args.list_exposures:
            print_summary(res)
            return

        # Select exposure (1-based index)
        try:
            exp = res.exposures[args.exposure - 1]
        except IndexError:
            sys.exit(f'[view_sca] Exposure index {args.exposure} out of range '
                     f'(found {len(res.exposures)} exposure(s); use --list to see them)')

        print(f'[view_sca] Selected exposure {args.exposure}/{len(res.exposures)}: '
              f'{exp.visit_id} exp {exp.exposure}', file=sys.stderr)

        # Validate SCA
        if args.sca is None:
            ap.error('--sca is required when using roman_mast query filters (--program, --pass, etc.)')

        if args.sca not in exp.scas:
            available = ', '.join(str(s) for s in sorted(exp.scas))
            sys.exit(f'[view_sca] SCA {args.sca} not in exposure. Available: {available}')

        sca_idx = exp.scas.index(args.sca)
        filename = exp.filenames[sca_idx]
        print(f'[view_sca] Streaming {filename}', file=sys.stderr)

        # Stream the single SCA
        asdf_file = res.missions.read_product(filename)
        try:
            data, detector, wcs_hdr, dq = stream_sca(
                None, None, asdf_file=asdf_file, sip_degree=args.sip_degree
            )
            title = (f'{detector} {exp.visit_id}  Exp: {exp.exposure:04d}  '
                     f'{exp.optical_element or "?"}')
        finally:
            if hasattr(asdf_file, 'close'):
                asdf_file.close()
    else:
        # Direct S3/local mode
        if args.uri is None or args.filename is None:
            ap.error('Either provide <uri> <filename> for direct S3/local access, '
                     'or use --program/--pass/etc. for roman_mast query mode')

        data, detector, wcs_hdr, dq = stream_sca(
            args.uri, args.filename, sip_degree=args.sip_degree
        )
        meta = parse_filename(args.filename)
        title = (f'{detector} {meta["visit_id"]}  '
                 f'Exp: {meta["exposure_num"]} '
                 f'{meta["filter"]}')

    dq_wanted = not args.no_dq

    # Delay photometry imports until actually needed to avoid hanging at module load
    from photometry import fit_background, run_aperture_photometry

    # Background subtraction (standalone — when --bkg without --phot)
    bkg_fit = data_sub = bkg_level = bkg_rms = None
    if args.bkg and not args.phot:
        bkg_fit, data_sub, _, bkg_level, bkg_rms = fit_background(
            data,
            dq=dq if dq_wanted else None,
            poly_degree=args.bkg_poly_degree,
        )

    # Run photometry if requested (includes background fitting internally)
    sources, phot_table, sp, stats = None, None, None, None
    if args.phot:
        sources, phot_table, sp, stats = run_aperture_photometry(
            data,
            dq=dq if dq_wanted else None,
            fwhm_pix=args.fwhm,
            detection_sigma=args.det_sigma,
            aperture_radius_fwhm=args.aper,
            annulus_inner_fwhm=args.annulus_inner,
            annulus_outer_fwhm=args.annulus_outer,
            bkg_poly_degree=args.bkg_poly_degree,
            snr_threshold=args.snr_threshold,
        )
        if phot_table is not None and args.phot_out:
            phot_table.write(args.phot_out, format='ascii.csv', overwrite=True)
            print(f'[view_sca] Photometry written to {args.phot_out}', file=sys.stderr)
        # Expose bkg results for residuals display when --bkg --phot
        if args.bkg and stats is not None:
            bkg_fit = stats.get('bkg_fit')
            data_sub = stats.get('data_sub')
            bkg_level = stats.get('bkg_level', 0.0)
            bkg_rms = stats.get('bkg_rms_median') or (
                float(np.median(stats['bkg_rms'])) if stats.get('bkg_rms') is not None else None
            )

    if args.display == 'mpl':
        display_in_mpl(
            data, detector,
            dq=dq if dq_wanted else None,
            title=title,
            wcs_header=wcs_hdr,
            show_channels=args.channels,
            show_gridlines=args.grid,
            sources=sources,
            phot_table=phot_table,
            sp=sp,
            stats=stats,
        )
        if args.bkg and data_sub is not None:
            display_residuals_mpl(
                data_sub, detector,
                bkg_fit=bkg_fit,
                bkg_level=bkg_level or 0.0,
                residual_rms=bkg_rms,
                poly_degree=args.bkg_poly_degree,
                dq=dq if dq_wanted else None,
                wcs_header=wcs_hdr,
                n_sigma=args.bkg_scale,
            )
    else:
        regions = make_channel_regions(
            data, detector,
            show_channels=args.channels,
            show_gridlines=args.grid,
        )
        display_in_ds9(
            data, regions,
            title=title,
            wcs_header=wcs_hdr,
            dq=dq if dq_wanted else None,
            ds9_target=args.ds9,
            scale=args.scale,
            cmap=args.cmap,
            sources=sources,
            phot_table=phot_table,
            sp=sp,
            data_sub=data_sub if args.bkg else None,
        )


if __name__ == '__main__':
    main()
