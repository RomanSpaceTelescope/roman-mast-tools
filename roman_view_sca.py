#!/usr/bin/env python
"""Stream a single Roman WFI SCA from S3 into DS9 with H4RG channel grid overlay.

DS9 must already be running before you call this script.

Usage
-----
    conda run -n roman-mast-tools python view_sca_ds9.py <uri> <filename>

    # public tutorial data (anonymous S3):
    python roman_view_sca.py \\
        s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \\
        r0003201001001001004_0001_wfi11_f106_cal.asdf

    # suppress channel dividers, keep section grid:
    python roman_view_sca.py --no-channels <uri> <filename>

    # connect to a named DS9 instance:
    python roman_view_sca.py --ds9 myds9 <uri> <filename>
"""

import argparse
import io
import sys

import numpy as np
import s3fs
import roman_datamodels as rdm
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

def stream_sca(uri, filename, *, sip_degree=4):
    """Open an ASDF file from S3 (anonymous) and return (data_f32, detector, wcs_hdr).

    Materialises ``dm.data`` and computes the SIP WCS approximation while the
    ASDF file handle is still open (gwcs needs the live tree for to_fits_sip).
    """
    path = uri.rstrip('/') + '/' + filename
    fs = s3fs.S3FileSystem(anon=True)
    print(f'[view_sca_ds9] streaming {path}', file=sys.stderr)
    with fs.open(path, 'rb') as f:
        dm = rdm.open(f)
        data     = np.asarray(dm.data[...], dtype=np.float32)
        detector = str(dm.meta.instrument.detector).upper()
        try:
            dq = np.asarray(dm.dq[...], dtype=np.int32)
        except AttributeError:
            dq = None   # uncal products have no DQ layer
        print(f'[view_sca_ds9] computing SIP WCS (degree={sip_degree}) ...', file=sys.stderr)
        gwcs    = dm.meta.wcs
        wcs_hdr = gwcs.to_fits_sip(bounding_box=gwcs.bounding_box, degree=sip_degree)
    print(f'[view_sca_ds9] {detector}  shape={data.shape}', file=sys.stderr)
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


# ---------------------------------------------------------------------------
# DS9 display
# ---------------------------------------------------------------------------

def display_in_ds9(data, regions, *, title='Roman WFI', wcs_header=None,
                   dq=None, ds9_target=None, scale='zscale', cmap='viridis'):
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

    # Build header: start from the SIP WCS if provided, then add OBJECT
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

    print('[view_sca_ds9] done — DS9 is open', file=sys.stderr)
    return d


# ---------------------------------------------------------------------------
# Matplotlib display
# ---------------------------------------------------------------------------

def display_in_mpl(data, detector, *, dq=None, title=None, wcs_header=None,
                   show_channels=False, show_gridlines=False,
                   figsize=(9, 8), norm=None):
    """Display a single SCA in an interactive matplotlib window.

    When wcs_header is provided the axes use WCSAxes projection so ticks and
    the status bar show RA/Dec; the data value (and DQ flag if present) are
    appended to the WCSAxes format_coord output.
    """
    import matplotlib.pyplot as plt
    from astropy.visualization import simple_norm

    ny, nx  = data.shape
    ref_pix = (4096 - nx) // 2
    sca_num = int(detector[-2:])
    rotated = _SCA_ROTATION[sca_num] == 180

    if norm is None:
        norm = simple_norm(data, 'asinh', vmin=0.5, vmax=4)

    projection = WCS(wcs_header) if wcs_header is not None else None
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

    _wcs_fmt = ax.format_coord   # WCSAxes provides RA/Dec; plain axes gives x/y
    def _fmt(x, y):
        col, row = int(round(x)), int(round(y))
        val_str = ''
        if 0 <= row < ny and 0 <= col < nx:
            val = data[row, col]
            dq_str = f'  dq={dq[row, col]}' if dq is not None else ''
            val_str = f'  x={col + ref_pix}  y={row + ref_pix}  val={val:.4g}{dq_str}'
        return _wcs_fmt(x, y) + val_str
    ax.format_coord = _fmt

    plt.tight_layout()
    plt.show()
    return fig, ax


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=(
            'Stream a single Roman WFI SCA from S3 and display it '
            'in DS9 or matplotlib with H4RG channel/section grid overlay.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('uri',      help='S3 base URI ending in /  (anonymous access)')
    ap.add_argument('filename', help='ASDF filename, e.g. r0003..._wfi11_f106_cal.asdf')
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
    args = ap.parse_args()

    dq_wanted = not args.no_dq
    data, detector, wcs_hdr, dq = stream_sca(args.uri, args.filename,
                                               sip_degree=args.sip_degree)

    meta = parse_filename(args.filename)
    title = (f'Roman WFI {detector}  Visit: {meta["visit_id"]}  '
             f'Exposure: {meta["exposure_num"]}  SCA: {meta["sca_num"]:02d}  '
             f'Filter: {meta["filter"]}')

    if args.display == 'mpl':
        display_in_mpl(
            data, detector,
            dq=dq if dq_wanted else None,
            title=title,
            wcs_header=wcs_hdr,
            show_channels=args.channels,
            show_gridlines=args.grid,
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
        )


if __name__ == '__main__':
    main()
