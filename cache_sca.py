#!/usr/bin/env python
"""Download and cache a single Roman WFI SCA ASDF file for local development.

Usage
-----
    python cache_sca.py <s3_uri> <filename> [--cache-dir ./cache]

Example
-------
    # Download from public tutorial data
    python cache_sca.py \\
        s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \\
        r0003201001001001004_0001_wfi11_f106_cal.asdf

    # Use a custom cache directory
    python cache_sca.py \\
        s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \\
        r0003201001001001004_0001_wfi11_f106_cal.asdf \\
        --cache-dir ~/roman_cache

Once cached, use with roman_view_sca.py:
    python roman_view_sca.py ./cache/ r0003201001001001004_0001_wfi11_f106_cal.asdf --display mpl
"""

import argparse
import os
import sys
import shutil

import s3fs


def main():
    ap = argparse.ArgumentParser(
        description='Download and cache a Roman WFI SCA ASDF file for local development.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('uri', help='S3 base URI ending in /')
    ap.add_argument('filename', help='ASDF filename')
    ap.add_argument('--cache-dir', default='./cache', metavar='PATH',
                    help='Local cache directory (default: ./cache)')
    args = ap.parse_args()

    cache_dir = os.path.expanduser(args.cache_dir)
    os.makedirs(cache_dir, exist_ok=True)

    s3_path = args.uri.rstrip('/') + '/' + args.filename
    local_path = os.path.join(cache_dir, args.filename)

    print(f'Downloading from S3: {s3_path}', file=sys.stderr)
    print(f'Cache destination:  {local_path}', file=sys.stderr)

    fs = s3fs.S3FileSystem(anon=True)
    try:
        fs.get(s3_path, local_path)
    except Exception as e:
        sys.exit(f'ERROR: Failed to download {s3_path}: {e}')

    file_size = os.path.getsize(local_path) / 1e6
    print(f'[cache_sca] Cached {args.filename} ({file_size:.1f} MB)', file=sys.stderr)
    print(f'\nTo use with roman_view_sca.py:\n'
          f'  python roman_view_sca.py {cache_dir} {args.filename} --display mpl\n',
          file=sys.stderr)


if __name__ == '__main__':
    main()
