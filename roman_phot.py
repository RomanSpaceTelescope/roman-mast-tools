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
"""

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from astropy.table import vstack

from roman_view_sca import stream_sca, run_aperture_photometry, parse_filename


# ---------------------------------------------------------------------------
# Per-SCA photometry
# ---------------------------------------------------------------------------

def phot_one_sca(uri, filename, *, phot_kwargs):
    """Stream one SCA and run aperture photometry.

    Returns a dict with keys: sca, detector, table, stats, sp.
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
        f'bkg={stats["bkg_level"]:.3g}  rms={stats["bkg_rms"]:.3g}  '
        f'n_sources={stats["n_sources"]}' +
        (f'  snr=[{float(snr_arr.min()):.1f}, {float(np.median(snr_arr)):.1f}, {float(snr_arr.max()):.1f}]'
         if snr_arr is not None and len(snr_arr) > 0 else ''),
        file=sys.stderr,
    )

    return {
        'sca': sca_num,
        'detector': detector,
        'table': phot_table,
        'stats': stats,
        'sp': sp,
    }


# ---------------------------------------------------------------------------
# Parallel fan-out
# ---------------------------------------------------------------------------

def phot_exposure(uri_filename_pairs, *, phot_kwargs=None, max_workers=8):
    """Run photometry on a list of (uri, filename) pairs in parallel.

    Returns a list of result dicts (sorted by SCA number), with None entries
    removed.
    """
    if phot_kwargs is None:
        phot_kwargs = {}

    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(phot_one_sca, uri, fn, phot_kwargs=phot_kwargs): fn
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

    results = phot_exposure(pairs, phot_kwargs=phot_kwargs, max_workers=args.workers)

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


if __name__ == '__main__':
    main()
