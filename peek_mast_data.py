"""Peek at available Roman WFI data on MAST.

Mirrors the search idiom from comm_streaming_example.ipynb:

    search = {'program': 114, 'pass': 57, 'detector': 'WFI04'}
    results  = missions.query_criteria(**search, select_cols=col_list)
    products = missions.get_unique_product_list(results)
    filtered = missions.filter_products(products, file_suffix='_cal')

The CLI exposes those same knobs as flags (--program, --pass, --detector,
--visit-id, --exposure) and lets you choose the file suffix via --data-level
(1 → '_uncal', 2 → '_cal', 'gw' → '_gw'). All criteria are optional; blank
values are simply omitted from the search dict.
"""

import keyring, keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

import warnings
warnings.filterwarnings('ignore')

import argparse

from streaming_utils import create_missions, get_MAST_token


# Columns to request from MAST — matches the notebook example
COL_LIST = [
    'program', 'execution_plan', 'pass', 'segment', 'visit', 'observation',
    'optical_element', 'exposure_type', 'instrument_name', 'detector',
    'productLevel', 'product_type', 'exposure_time',
    'exposure_start_time', 'exposure_end_time', 'fileSetName',
]

# Map data-level shorthand → filter_products file_suffix
LEVEL_SUFFIX = {1: '_uncal', 2: '_cal', 'gw': '_gw'}


def build_search(*, program=None, pass_=None, detector=None,
                 visit_id=None, exposure=None):
    """Assemble a MAST search dict from the given criteria, omitting blanks.

    Any criterion left as None (or empty string) is not sent to MAST — so
    calling with no arguments searches everything, exactly like the notebook
    when you leave keys out of ``search``.
    """
    search = {}
    if program not in (None, ''):
        search['program'] = int(program)
    if pass_ not in (None, ''):
        search['pass'] = int(pass_)
    if detector not in (None, ''):
        search['detector'] = str(detector).upper()
    if visit_id not in (None, ''):
        search['visit_id'] = visit_id
    if exposure not in (None, ''):
        # observation_id ends with the 4-digit exposure number; MAST accepts
        # '*NNNN' as a wildcard suffix match (see streaming_utils.group_file_urls).
        search['observation_id'] = f'*{int(exposure):04d}'
    return search


def search_products(missions, *, search, data_level=2):
    """Run the notebook query → unique products → filter_products pipeline.

    Returns the ``filtered`` astropy Table (may be empty).
    """
    print(f"Search criteria: {search or '(none — will return all data)'}")
    results = missions.query_criteria(**search, select_cols=COL_LIST)
    print(f"Total number of results: {len(results)}")
    if len(results) == 0:
        return results

    products = missions.get_unique_product_list(results)

    suffix = LEVEL_SUFFIX.get(data_level)
    if suffix is None:
        raise ValueError(f"Unknown data_level {data_level!r}; expected one of {list(LEVEL_SUFFIX)}")

    filtered = missions.filter_products(products, file_suffix=suffix)
    print(f"After filter_products(file_suffix={suffix!r}): {len(filtered)} products")
    return filtered


def print_product_table(filtered, max_rows=50):
    """Print the filenames in the filtered product table, one per line."""
    if len(filtered) == 0:
        print("\n(no products match)")
        return

    print(f"\nProducts ({len(filtered)}):")
    for i, name in enumerate(filtered['filename']):
        if i >= max_rows:
            print(f"  ... and {len(filtered) - max_rows} more")
            break
        print(f"  {name}")


def peek_mast_data(*, program=None, pass_=None, detector=None,
                   visit_id=None, exposure=None, data_level=2,
                   mast_token=None, max_rows=50):
    """Convenience wrapper: authenticate, search, filter, print.

    Returns the filtered astropy Table so callers can drive downstream code.
    """
    if mast_token is None:
        mast_token = get_MAST_token()
    if not mast_token:
        raise ValueError("MAST_API_TOKEN not found in .env file")

    missions = create_missions(mast_token)

    search = build_search(program=program, pass_=pass_, detector=detector,
                          visit_id=visit_id, exposure=exposure)
    filtered = search_products(missions, search=search, data_level=data_level)
    print_product_table(filtered, max_rows=max_rows)
    return filtered


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _prompt(msg, default=None):
    hint = f" (default {default!r})" if default not in (None, '') else ''
    val = input(f"{msg}{hint}: ").strip()
    return val if val else default


def interactive():
    """Prompt the user for each search criterion (all optional)."""
    print("\nEnter search criteria (blank = don't filter on this field).")
    print("Example from the notebook: program=114, pass=57, detector=WFI04\n")

    program   = _prompt("Program (int)")
    pass_     = _prompt("Pass (int)")
    detector  = _prompt("Detector (e.g. WFI04)")
    visit_id  = _prompt("Visit ID (e.g. 0011401057001001001)")
    exposure  = _prompt("Exposure number (int)")

    lvl_raw = _prompt("Data level (1=uncal, 2=cal, gw=guide-window)", '2')
    data_level = 'gw' if str(lvl_raw).lower() == 'gw' else int(lvl_raw)

    return dict(program=program, pass_=pass_, detector=detector,
                visit_id=visit_id, exposure=exposure, data_level=data_level)


def parse_level(s):
    if s is None:
        return 2
    if str(s).lower() == 'gw':
        return 'gw'
    return int(s)


def main():
    parser = argparse.ArgumentParser(
        description='Peek at available Roman WFI data on MAST (notebook-style search).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python peek_mast_data.py --program 114 --pass 57 --detector WFI04
  python peek_mast_data.py --visit-id 0011401057001001001 --exposure 4
  python peek_mast_data.py --program 114 --data-level 1
  python peek_mast_data.py                       (interactive)
        """,
    )
    parser.add_argument('--program',   type=int, default=None, help='APT program ID (e.g. 114)')
    parser.add_argument('--pass',      dest='pass_', type=int, default=None,
                        metavar='PASS', help='Pass number (e.g. 57)')
    parser.add_argument('--detector',  default=None, help='Detector, e.g. WFI04')
    parser.add_argument('--visit-id',  default=None, help='Visit ID, e.g. 0011401057001001001')
    parser.add_argument('--exposure',  type=int, default=None,
                        help='Exposure number (matches trailing 4 digits of observation_id)')
    parser.add_argument('--data-level', default='2',
                        help="1=uncal, 2=cal, gw=guide-window (default 2)")
    parser.add_argument('--max-rows',  type=int, default=50,
                        help='Cap on rows printed (default 50)')

    args = parser.parse_args()

    any_flag = any(v is not None for v in [
        args.program, args.pass_, args.detector, args.visit_id, args.exposure,
    ])

    try:
        if not any_flag:
            print("Roman WFI MAST Data Peek")
            print("=" * 70)
            params = interactive()
            peek_mast_data(max_rows=args.max_rows, **params)
        else:
            peek_mast_data(
                program=args.program,
                pass_=args.pass_,
                detector=args.detector,
                visit_id=args.visit_id,
                exposure=args.exposure,
                data_level=parse_level(args.data_level),
                max_rows=args.max_rows,
            )
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
