"""Peek at available data on MAST based on visit/exposure/SCA wildcard filtering.

Shows what data is available without actually downloading or extracting metadata.
"""

import keyring, keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

import warnings
warnings.filterwarnings('ignore')

import argparse
from query_utils import (
    add_query_args,
    prompt_query_params,
    resolve_query,
)


def peek_mast_data(visit_wildcard='*', exp_wildcard='*', sca_wildcard='*', data_level=2, verbose=False):
    """Peek at available data on MAST based on wildcard filters.

    Args:
        visit_wildcard (str): Visit ID pattern with wildcards (e.g., '001240100100100*')
        exp_wildcard (str): Exposure specifier (e.g., '1', '1-5', '*')
        sca_wildcard (str): SCA specifier (e.g., '1', '01,02,03', '01-10', '*')
        data_level (int): Data level (1=uncal, 2=cal)
        verbose (bool): Show detailed output (currently unused, reserved for future use)

    Returns:
        dict: Summary of available data keyed by visit ID
    """
    q = resolve_query(
        visit_id=visit_wildcard,
        exp_spec=exp_wildcard if exp_wildcard != '*' else None,
        data_level=data_level,
        sca_spec=sca_wildcard if sca_wildcard not in ('*', 'all') else None,
    )

    # Build the legacy return format expected by callers
    summary = {
        'total_scas': len(q.scas),
        'filtered_scas': len(q.scas),
        'data_by_visit': {
            vid: {'exposures': info['exposures'], 'scas': info['scas']}
            for vid, info in q.summary.items()
        },
        'total_files': sum(
            len(info['exposures']) * len(info['scas'])
            for info in q.summary.values()
        ),
    }

    print(f"You can now use export_metadata_csv.py with these parameters:")
    print(f"  visit_id='{visit_wildcard}'")
    if q.exp_nums:
        print(f"  exp_nums={q.exp_nums}")
    print(f"  data_level={data_level}")

    return summary


def main():
    """Interactive or command-line entry point."""
    parser = argparse.ArgumentParser(
        description='Peek at available Roman WFI data on MAST.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python peek_mast_data.py '*'
  python peek_mast_data.py '00124010*' --exp-range '1-3'
  python peek_mast_data.py '0012401001001002001' --data-level 1 --scas '1,2,3'
  python peek_mast_data.py  (interactive mode)
        """,
    )
    add_query_args(parser, visit_wildcard=True, exp_mode='flexible', sca_mode='all')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()

    if args.visit_id is None:
        print("MAST Data Peek")
        print("=" * 70)
        params = prompt_query_params(
            visit_wildcard=True,
            exp_mode='flexible',
            sca_mode='all',
            defaults={'visit_id': '*'},
        )
        visit_id = params['visit_id']
        exp_spec = params['exp_spec']
        data_level = params['data_level']
        sca_spec = params['sca_spec']
        verbose = False
    else:
        visit_id = args.visit_id
        exp_spec = getattr(args, 'exp_range', None) or (str(args.exp_num) if getattr(args, 'exp_num', None) else None)
        data_level = args.data_level if args.data_level is not None else 2
        sca_spec = args.scas
        verbose = args.verbose

    try:
        peek_mast_data(
            visit_wildcard=visit_id,
            exp_wildcard=exp_spec or '*',
            sca_wildcard=sca_spec or '*',
            data_level=data_level,
            verbose=verbose,
        )
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
