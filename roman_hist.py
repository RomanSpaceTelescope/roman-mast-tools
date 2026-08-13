#!/usr/bin/env python
"""Make histograms from roman_phot.py CSV output.

Reads a CSV file from aperture photometry and generates histograms
of selected columns, focusing on aperture_sum_0 and SNR by default.

Usage
-----
    # Basic usage with default columns (aperture_sum_0 and snr):
    python roman_hist.py roman_phot_combined.csv

    # Custom columns and output:
    python roman_hist.py roman_phot_combined.csv \\
        --columns aperture_sum_0 flux_bkgsub \\
        --output histograms.png \\
        --bins 50

    # Generate histograms for all numeric columns:
    python roman_hist.py roman_phot_combined.csv \\
        --all-columns --output all_hists.png

    # Aggressive outlier rejection (IQR × 1.0):
    python roman_hist.py roman_phot_combined.csv \\
        --outlier-method iqr --iqr-mult 1.0 --output hists_clean.png

    # Sigma clipping (keep ±2σ):
    python roman_hist.py roman_phot_combined.csv \\
        --outlier-method sigma --sigma-clip 2.0 --output hists_sigma.png

    # Disable outlier rejection:
    python roman_hist.py roman_phot_combined.csv \\
        --outlier-method none --output hists_raw.png
"""

import argparse
import sys

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table


def load_csv(path):
    """Load CSV file using astropy Table."""
    try:
        table = Table.read(path, format='ascii.csv')
        return table
    except Exception as exc:
        print(f'[roman_hist] ERROR: failed to read {path}: {exc}', file=sys.stderr)
        sys.exit(1)


def get_numeric_columns(table, exclude_spatial=False):
    """Get all numeric column names from table."""
    numeric = []
    for colname in table.colnames:
        col = table[colname]
        try:
            # Check if column is numeric
            if np.issubdtype(col.dtype, np.number):
                # Optionally exclude spatial/ID columns
                if exclude_spatial and colname in ('sca', 'detector', 'id', 'x_center', 'y_center'):
                    continue
                numeric.append(colname)
        except:
            pass
    return numeric


def validate_columns(table, columns, all_columns=False):
    """Validate requested columns exist in table."""
    if all_columns:
        cols = get_numeric_columns(table, exclude_spatial=True)
        if not cols:
            print('[roman_hist] ERROR: no numeric columns found in table', file=sys.stderr)
            sys.exit(1)
        return cols

    # Check requested columns exist
    missing = [c for c in columns if c not in table.colnames]
    if missing:
        print(f'[roman_hist] ERROR: columns not found: {missing}', file=sys.stderr)
        print(f'[roman_hist] Available columns: {table.colnames}', file=sys.stderr)
        sys.exit(1)

    # Check columns are numeric
    non_numeric = []
    for c in columns:
        col = table[c]
        if not np.issubdtype(col.dtype, np.number):
            non_numeric.append(c)
    if non_numeric:
        print(f'[roman_hist] ERROR: non-numeric columns: {non_numeric}', file=sys.stderr)
        sys.exit(1)

    return columns


def reject_outliers(data, method='iqr', iqr_mult=1.5, sigma_clip=3):
    """Reject outliers from data using specified method.

    Parameters
    ----------
    data : array-like
        Input data
    method : {'iqr', 'sigma', 'none'}
        Rejection method:
        - 'iqr': Interquartile range (robust to outliers)
        - 'sigma': Sigma clipping (assumes normal distribution)
        - 'none': No rejection
    iqr_mult : float
        IQR multiplier for outlier threshold (default 1.5 = standard Tukey method)
    sigma_clip : float
        Number of standard deviations for sigma clipping

    Returns
    -------
    filtered_data : array
        Data with outliers removed
    n_rejected : int
        Number of points rejected
    """
    if method == 'none':
        return data, 0

    if method == 'iqr':
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        lower = q1 - iqr_mult * iqr
        upper = q3 + iqr_mult * iqr
        mask = (data >= lower) & (data <= upper)
    elif method == 'sigma':
        mean = np.mean(data)
        std = np.std(data)
        mask = np.abs(data - mean) <= sigma_clip * std
    else:
        raise ValueError(f'Unknown outlier method: {method}')

    n_rejected = len(data) - np.sum(mask)
    return data[mask], n_rejected


def make_histograms(table, columns, bins=30, figsize=None, title=None,
                   outlier_method='iqr', outlier_iqr_mult=1.5, outlier_sigma=3):
    """Create histogram figure for specified columns.

    Parameters
    ----------
    table : astropy.table.Table
        Input data table
    columns : list of str
        Column names to histogram
    bins : int
        Number of bins per histogram
    figsize : tuple, optional
        Figure size (width, height)
    title : str, optional
        Overall figure title
    outlier_method : {'iqr', 'sigma', 'none'}
        Outlier rejection method
    outlier_iqr_mult : float
        IQR multiplier for outlier rejection
    outlier_sigma : float
        Sigma clip threshold

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    n_cols = len(columns)
    if n_cols == 0:
        raise ValueError('No columns to plot')

    # Layout: 2 columns, as many rows as needed
    n_rows = (n_cols + 1) // 2
    if figsize is None:
        figsize = (12, 4 * n_rows)

    fig, axes = plt.subplots(n_rows, 2, figsize=figsize)
    if n_cols == 1:
        axes = np.array([[axes, axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)

    axes_flat = axes.flatten()

    for idx, colname in enumerate(columns):
        ax = axes_flat[idx]
        data = table[colname]

        # Filter out non-finite values
        valid_data = data[np.isfinite(data)]
        if len(valid_data) == 0:
            ax.text(0.5, 0.5, f'No valid data in {colname}',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(colname)
            continue

        # Reject outliers
        plot_data = valid_data
        n_rejected = 0
        if outlier_method != 'none':
            plot_data, n_rejected = reject_outliers(
                valid_data, method=outlier_method,
                iqr_mult=outlier_iqr_mult, sigma_clip=outlier_sigma
            )

        if len(plot_data) == 0:
            ax.text(0.5, 0.5, f'No data left after outlier rejection in {colname}',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(colname)
            continue

        ax.hist(plot_data, bins=bins, edgecolor='black', alpha=0.7)
        ax.set_xlabel(colname)
        ax.set_ylabel('Count')

        title_str = f'{colname} (n={len(plot_data)}'
        if n_rejected > 0:
            title_str += f', {n_rejected} outliers removed'
        title_str += ')'
        ax.set_title(title_str)
        ax.grid(True, alpha=0.3)

        # Add statistics text
        stats_text = f'μ={np.mean(plot_data):.3g}\nσ={np.std(plot_data):.3g}\nmed={np.median(plot_data):.3g}'
        ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Hide extra subplots if odd number of columns
    for idx in range(len(columns), len(axes_flat)):
        axes_flat[idx].axis('off')

    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold')
    fig.tight_layout()

    return fig


def main():
    ap = argparse.ArgumentParser(
        description='Make histograms from roman_phot.py CSV output.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        'csv_file', metavar='CSV_FILE',
        help='Input CSV file from roman_phot.py',
    )
    ap.add_argument(
        '-c', '--columns', nargs='+', default=['aperture_sum_0', 'snr'],
        metavar='COL',
        help='Columns to histogram (default: aperture_sum_0 snr)',
    )
    ap.add_argument(
        '-a', '--all-columns', action='store_true',
        help='Histogram all numeric columns (except spatial/ID)',
    )
    ap.add_argument(
        '-b', '--bins', type=int, default=30, metavar='N',
        help='Number of bins per histogram (default: 30)',
    )
    ap.add_argument(
        '-o', '--output', metavar='PATH',
        help='Output image file (e.g., histograms.png, histograms.pdf)',
    )
    ap.add_argument(
        '-t', '--title', metavar='TITLE',
        help='Overall figure title',
    )
    ap.add_argument(
        '--width', type=float, default=12, metavar='N',
        help='Figure width in inches (default: 12)',
    )
    ap.add_argument(
        '--height-per-row', type=float, default=4, metavar='N',
        help='Figure height per row in inches (default: 4)',
    )
    ap.add_argument(
        '--outlier-method', choices=['iqr', 'sigma', 'none'], default='iqr',
        metavar='METHOD',
        help='Outlier rejection method: iqr (Tukey fences, default), sigma (σ-clipping), or none',
    )
    ap.add_argument(
        '--iqr-mult', type=float, default=1.5, metavar='N',
        help='IQR multiplier for outlier threshold (default: 1.5 = standard Tukey method)',
    )
    ap.add_argument(
        '--sigma-clip', type=float, default=3.0, metavar='N',
        help='Number of σ for sigma clipping (default: 3.0)',
    )

    args = ap.parse_args()

    # Load CSV
    print(f'[roman_hist] Reading {args.csv_file}...', file=sys.stderr)
    table = load_csv(args.csv_file)
    print(f'[roman_hist] Loaded {len(table)} rows, {len(table.colnames)} columns',
          file=sys.stderr)

    # Validate and get columns
    columns = validate_columns(table, args.columns, all_columns=args.all_columns)
    print(f'[roman_hist] Plotting: {columns}', file=sys.stderr)

    # Calculate figure size
    n_rows = (len(columns) + 1) // 2
    figsize = (args.width, args.height_per_row * n_rows)

    # Create histograms
    title = args.title
    if not title and not args.all_columns:
        title = f'Histograms: {", ".join(columns)}'

    fig = make_histograms(
        table, columns, bins=args.bins, figsize=figsize, title=title,
        outlier_method=args.outlier_method,
        outlier_iqr_mult=args.iqr_mult,
        outlier_sigma=args.sigma_clip,
    )

    # Save or show
    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches='tight')
        print(f'[roman_hist] Saved to {args.output}', file=sys.stderr)
    else:
        plt.show()


if __name__ == '__main__':
    main()
