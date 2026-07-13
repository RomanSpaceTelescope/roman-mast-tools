"""Shared query-parameter helpers for roman-mast-tools scripts.

Provides:
  parse_exp_spec / parse_sca_spec  — parse user strings into int lists
  add_query_args                   — attach standard argparse args to a parser
  prompt_query_params              — interactive prompts returning a param dict
  ResolvedQuery / resolve_query    — MAST lookup + filtering + summary display
"""

import re
import argparse
from dataclasses import dataclass, field

from streaming_utils import get_MAST_token, group_file_urls, create_missions, get_exp_num, get_sca_num


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_exp_spec(s):
    """Parse an exposure specifier into a sorted list of ints, or None (=all).

    Accepts:
      None / '' / '*'   → None (caller uses all available)
      '1'               → [1]
      '1-3'             → [1, 2, 3]
      '1,2,5'           → [1, 2, 5]
      '1-3,5,7-8'       → [1, 2, 3, 5, 7, 8]
    """
    if s is None or str(s).strip() in ('', '*'):
        return None

    s = str(s).strip()
    nums = set()

    for part in s.split(','):
        part = part.strip()
        if '-' in part:
            lo, hi = part.split('-', 1)
            nums.update(range(int(lo.strip()), int(hi.strip()) + 1))
        else:
            nums.add(int(part))

    return sorted(nums)


def parse_sca_spec(s):
    """Parse an SCA specifier into a sorted list of ints, or None (=all).

    Accepts the same syntax as parse_exp_spec, plus 'all'.
    """
    if s is None or str(s).strip().lower() in ('', '*', 'all'):
        return None

    return parse_exp_spec(s)


# ---------------------------------------------------------------------------
# Argparse integration
# ---------------------------------------------------------------------------

def add_query_args(parser, *, visit_wildcard=False, exp_mode='single', sca_mode='all'):
    """Add standard Roman query arguments to *parser*.

    Parameters
    ----------
    visit_wildcard : bool
        If True, label the positional as a wildcard pattern (e.g. peek_mast_data).
    exp_mode : str
        'single'   — --exp-num only (int)
        'range'    — --exp-range only (string)
        'flexible' — --exp-num + --exp-range (mutually exclusive group)
    sca_mode : str
        'single'   — --sca only (single int)
        'multi'    — --scas (comma/range string)
        'all'      — --scas (comma/range string, default all)
    """
    visit_help = (
        'Visit ID wildcard pattern (e.g. "001240100100100*")' if visit_wildcard
        else 'Roman visit ID (e.g. "0012401001001002001")'
    )
    parser.add_argument('visit_id', nargs='?', default=None, help=visit_help)

    parser.add_argument(
        '--data-level', type=int, default=None, choices=[1, 2],
        help='Data level: 1=uncal, 2=cal (default: 2)',
    )

    if exp_mode == 'single':
        parser.add_argument('--exp-num', type=int, default=None,
                            help='Exposure number within visit (default: 1)')
    elif exp_mode == 'range':
        parser.add_argument('--exp-range', default=None,
                            help='Exposure range, e.g. "1-3" or "1,2,5" (default: all)')
    elif exp_mode == 'flexible':
        grp = parser.add_mutually_exclusive_group()
        grp.add_argument('--exp-num', type=int, default=None,
                         help='Single exposure number')
        grp.add_argument('--exp-range', default=None,
                         help='Exposure range, e.g. "1-3" or "1,2,5" (default: all)')

    if sca_mode == 'single':
        parser.add_argument('--sca', type=int, default=None,
                            help='SCA number 1-18 (default: 1)')
    else:
        parser.add_argument('--scas', default=None,
                            help='SCA numbers, e.g. "1,2,3" or "1-5" (default: all)')


# ---------------------------------------------------------------------------
# Interactive prompting
# ---------------------------------------------------------------------------

def prompt_query_params(
    *,
    visit_wildcard=False,
    exp_mode='flexible',
    sca_mode='all',
    defaults=None,
):
    """Interactively prompt the user for query parameters.

    Returns a dict with keys:
      visit_id   (str)
      exp_spec   (str — raw, pass to parse_exp_spec / resolve_query)
      data_level (int)
      sca_spec   (str — raw, pass to parse_sca_spec / resolve_query)
    """
    if defaults is None:
        defaults = {}

    # --- Visit ---
    if visit_wildcard:
        print("\nVisit wildcard examples:")
        print("  '*'                  — all data")
        print("  '00124010*'          — all visits in program 124010")
        print("  '001240100100100*'   — narrow visit range")
        print("  '0012401001001002001'— exact visit")
        default_visit = defaults.get('visit_id', '*')
        val = input(f"\nEnter visit wildcard pattern (default {default_visit!r}): ").strip()
        visit_id = val if val else default_visit
    else:
        print("\nVisit ID examples:")
        print("  '0012401001001002001' — exact visit")
        print("  '001240100100100*'    — wildcard (matches multiple visits)")
        default_visit = defaults.get('visit_id', '0012401001001002001')
        val = input(f"\nEnter Visit ID (default {default_visit!r}): ").strip()
        visit_id = val if val else default_visit

    # --- Exposure ---
    if exp_mode == 'single':
        default_exp = str(defaults.get('exp_num', 1))
        val = input(f"Enter exposure number (1-based, default {default_exp}): ").strip()
        exp_spec = val if val else default_exp
    elif exp_mode == 'range':
        print("\nExposure examples:  '1-3'  '1,2,5'  (blank = all)")
        val = input("Enter exposure range (default all): ").strip()
        exp_spec = val if val else None
    else:  # flexible
        print("\nExposure examples:  '1'  '1-3'  '1,2,5'  (blank = all)")
        val = input("Enter exposure(s) (default all): ").strip()
        exp_spec = val if val else None

    # --- Data level ---
    default_dl = str(defaults.get('data_level', 2))
    val = input(f"Enter data level (1=uncal, 2=cal, default {default_dl}): ").strip()
    data_level = int(val) if val else int(default_dl)

    # --- SCA ---
    if sca_mode == 'single':
        default_sca = str(defaults.get('sca', 1))
        val = input(f"Enter SCA number (1-18, default {default_sca}): ").strip()
        sca_spec = val if val else default_sca
    else:
        print("\nSCA examples:  '1'  '1,2,3'  '1-10'  (blank = all)")
        val = input("Enter SCA(s) (default all): ").strip()
        sca_spec = val if val else None

    return {
        'visit_id': visit_id,
        'exp_spec': exp_spec,
        'data_level': data_level,
        'sca_spec': sca_spec,
    }


# ---------------------------------------------------------------------------
# ResolvedQuery + resolve_query
# ---------------------------------------------------------------------------

@dataclass
class ResolvedQuery:
    """Result of a MAST query + filter operation, ready to pass to a core function."""
    visit_id: str
    exp_nums: list
    data_level: int
    scas: list
    urls: dict
    summary: dict = field(default_factory=dict)
    missions: object = None  # logged-in MastMissions session


def resolve_query(visit_id, exp_spec=None, data_level=2, sca_spec=None, mast_token=None):
    """Query MAST, apply filters, print a summary table, and return a ResolvedQuery.

    Parameters
    ----------
    visit_id   : str — visit ID or wildcard pattern
    exp_spec   : None / int / str — parsed by parse_exp_spec
    data_level : int — 1 or 2
    sca_spec   : None / int / str — parsed by parse_sca_spec
    mast_token : str or None — loaded from .env if None

    Returns
    -------
    ResolvedQuery with missions session for streaming operations
    """
    if mast_token is None:
        mast_token = get_MAST_token()
    if not mast_token:
        raise ValueError("MAST_API_TOKEN not found in .env file")

    print("Querying MAST for available data...")
    missions = create_missions(mast_token)
    missions_obj, raw_files = group_file_urls(missions=missions, visit_id=visit_id, exp_num='', data_level=data_level)

    if not raw_files:
        raise ValueError(f"No data found for visit pattern '{visit_id}' at level {data_level}")

    all_scas = sorted(raw_files.keys())

    # --- Parse filters ---
    exp_list = parse_exp_spec(exp_spec)
    sca_list = parse_sca_spec(sca_spec)

    scas = [s for s in all_scas if sca_list is None or s in sca_list]
    if not scas:
        raise ValueError(f"No SCAs match the requested filter (available: {all_scas})")

    # --- Build per-visit summary and collect filenames ---
    summary = {}
    matched_files = []
    available_exps = set()
    for sca_num in scas:
        for exp_num, filename in raw_files[sca_num].items():
            available_exps.add(exp_num)
            if exp_list is not None and exp_num not in exp_list:
                continue
            # Extract visit ID from filename (format: r{visit_id}_{exp}_{sca}_{filter}_{ext})
            m = re.search(r'r(\d{18,})', filename)
            if m:
                vid = m.group(1)
                if vid not in summary:
                    summary[vid] = {'exposures': set(), 'scas': set()}
                summary[vid]['exposures'].add(exp_num)
                summary[vid]['scas'].add(sca_num)
            matched_files.append(filename)

    if not summary:
        raise ValueError(
            f"No exposures match the requested filter {exp_spec!r} "
            f"(requested exp_list={exp_list}, available in SCAs {scas}: {sorted(available_exps)})"
        )

    exp_nums = sorted({e for info in summary.values() for e in info['exposures']})

    # --- Build {sca: {real_exp_num: filename}} for streaming callers ---
    filtered_files = {}
    for s in scas:
        filtered_files[s] = {}
        for exp_num, filename in raw_files[s].items():
            if exp_list is None or exp_num in exp_list:
                filtered_files[s][exp_num] = filename

    # --- Print summary table ---
    print("\n" + "=" * 70)
    print("AVAILABLE DATA")
    print("=" * 70)
    print(f"  Visit pattern : {visit_id}")
    print(f"  Data level    : {data_level}")
    print(f"  Exposures     : {exp_nums}")
    print(f"  SCAs          : {scas}")
    print()

    if summary:
        print(f"  {'Visit ID':<22} {'Exposures':>10} {'SCAs':>6}")
        print("  " + "-" * 40)
        for vid in sorted(summary):
            info = summary[vid]
            print(f"  {vid:<22} {len(info['exposures']):>10} {len(info['scas']):>6}")
        print("  " + "-" * 40)
    else:
        print("  (No visit IDs could be parsed from filenames)")

    if matched_files:
        print(f"\n  Files ({len(matched_files)}):")
        for fname in sorted(matched_files):
            print(f"    {fname}")

    print()

    return ResolvedQuery(
        visit_id=visit_id,
        exp_nums=exp_nums,
        data_level=data_level,
        scas=scas,
        urls=filtered_files,
        summary=summary,
        missions=missions_obj,
    )
