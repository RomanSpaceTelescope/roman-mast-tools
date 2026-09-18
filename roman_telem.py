#!/usr/bin/env python3
"""
roman_telem.py
--------------
Command-line and Python API for retrieving Roman Space Telescope telemetry
(mnemonics) from the MAST Engineering Database, with optional plotting via
`roman_telem_plot`.

Wraps `mast_eng_db_query.MASTEngDBQuery` and adds:
    * Positional group names as the shortest query form:
          roman-telem rcs_pd
          roman-telem bank1_leds bank2_leds --days 7
      Group names are looked up in tlm_groups.yaml (cwd) or $ROMAN_TELEM_GROUPS.
    * Mnemonic input via -m, --mnemonics-file, or --plot-groups YAML config
    * Group filtering (--select-groups) and discovery (--list-groups)
    * Time range via -s/-e or --days (default: last 1 day)
    * CSV/Parquet/HDF5/Pickle output via --output
    * Trending plots via --plot-output / --show, with --plot-per-mnemonic
      for one subplot per mnemonic
"""

from __future__ import annotations

import argparse
import io
import os
import sys
import time
import requests
from datetime import datetime
from typing import List, Optional, Sequence, Tuple, Union
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd


# ---------------------------------------------------------------------------
# MASTEngDBQuery (inlined from mast_eng_db_query.py)
# ---------------------------------------------------------------------------

class MASTEngDBQuery:
    def __init__(self, api_token: Optional[str] = None, server: str = 'mast'):
        if server == 'int':
            self.base_url = "https://mastint.stsci.edu/edp/api/v0.1/mnemonics/spa/roman/data"
            self.api_token = api_token or os.environ.get('MAST_API_TOKEN_INT')
        else:
            self.base_url = "https://mast.stsci.edu/edp/api/v0.1/mnemonics/spa/roman/data"
            self.api_token = api_token or os.environ.get('MAST_API_TOKEN')

        if not self.api_token:
            raise ValueError("API token must be provided or set in MAST_API_TOKEN environment variable")

    def _format_datetime(self, dt: Union[str, datetime]) -> str:
        if isinstance(dt, str):
            try:
                dt = pd.to_datetime(dt)
            except:
                return dt

        if isinstance(dt, datetime) or isinstance(dt, pd.Timestamp):
            return dt.strftime('%Y-%m-%dT%H:%M:%S')

        return str(dt)

    def query_mnemonic(
        self,
        mnemonic: str,
        start_time: Union[str, datetime],
        end_time: Union[str, datetime],
        result_format: str = "csv",
        verbose=True,
    ) -> pd.DataFrame:
        params = {
            "mnemonic": mnemonic,
            "s_time": self._format_datetime(start_time),
            "e_time": self._format_datetime(end_time),
            "result_format": result_format,
        }

        headers = {"Authorization": f"token {self.api_token}"}

        response = requests.get(self.base_url, params=params, headers=headers)
        response.raise_for_status()

        if result_format == "csv":
            df = pd.read_csv(io.StringIO(response.text))
        else:
            df = pd.DataFrame(response.json())

        df["mnemonic"] = mnemonic

        for col in df.columns:
            if "time" in col.lower() or "date" in col.lower():
                try:
                    df[col] = pd.to_datetime(df[col])
                except:
                    pass

        if verbose:
            print("Mnemonic " + mnemonic + " query completed.")
        return df

    def _query_mnemonic_with_retries(
        self,
        mnemonic: str,
        start_time: Union[str, datetime],
        end_time: Union[str, datetime],
        result_format: str,
        max_retries: int = 3,
        retry_delay: float = 1.0,
        verbose=True,
    ) -> pd.DataFrame:
        last_exception = None
        attempts = max(1, max_retries)

        for attempt in range(1, attempts + 1):
            try:
                return self.query_mnemonic(
                    mnemonic, start_time, end_time, result_format, verbose=verbose
                )
            except Exception as e:
                last_exception = e
                if attempt < attempts:
                    wait = retry_delay * attempt
                    if verbose:
                        print(
                            f"⚠ Attempt {attempt}/{attempts} failed for "
                            f"{mnemonic}: {str(e)}. Retrying in {wait:.1f}s..."
                        )
                    time.sleep(wait)
                else:
                    if verbose:
                        print(
                            f"✗ Error querying {mnemonic} after {attempts} "
                            f"attempts: {str(e)}"
                        )

        return pd.DataFrame()

    def query_multiple_mnemonics(
        self,
        mnemonics: List[str],
        start_time: Union[str, datetime],
        end_time: Union[str, datetime],
        combine: bool = True,
        result_format: str = "csv",
        pbarparent=None,
        max_retries: int = 3,
        retry_delay: float = 1.0,
        verbose=True
    ) -> Union[pd.DataFrame, dict]:
        results = {}

        if len(mnemonics) <= 1:
            mnemonic_count = 0
            for mnemonic in mnemonics:
                if pbarparent is not None:
                    pbarparent.progressbar_increment(
                        mnemonic_count, len(mnemonics), "Querying " + mnemonic
                    )
                results[mnemonic] = self._query_mnemonic_with_retries(
                    mnemonic,
                    start_time,
                    end_time,
                    result_format,
                    max_retries=max_retries,
                    retry_delay=retry_delay,
                    verbose=verbose,
                )
                mnemonic_count += 1
        else:
            if pbarparent is not None:
                pbarparent.progressbar_increment(
                    0, len(mnemonics), "Starting concurrent queries..."
                )

            max_workers = min(10, len(mnemonics))
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_mnemonic = {
                    executor.submit(
                        self._query_mnemonic_with_retries,
                        mnemonic,
                        start_time,
                        end_time,
                        result_format,
                        max_retries,
                        retry_delay,
                        verbose=verbose
                    ): mnemonic
                    for mnemonic in mnemonics
                }

                completed_count = 0
                for future in as_completed(future_to_mnemonic):
                    mnemonic = future_to_mnemonic[future]
                    try:
                        df = future.result()
                        results[mnemonic] = df
                    except Exception as e:
                        if verbose:
                            print(f"✗ Unexpected error for {mnemonic}: {str(e)}")
                        results[mnemonic] = pd.DataFrame()

                    completed_count += 1
                    if pbarparent is not None:
                        pbarparent.progressbar_increment(
                            completed_count,
                            len(mnemonics),
                            f"Completed {mnemonic}",
                        )

        if not combine:
            return results

        if pbarparent is not None:
            pbarparent.progressbar_increment(
                len(mnemonics), len(mnemonics), "Mnemonics Queried."
            )

        combined_dfs = []
        mnemonic_count = 0
        start_cpu = time.process_time()
        for mnemonic, df in results.items():
            if pbarparent is not None:
                pbarparent.progressbar_increment(
                    mnemonic_count,
                    len(mnemonics),
                    "Combine: Appending " + mnemonic,
                )
            if (
                not df.empty
                and "ObsTime" in df.columns
                and "EUValue" in df.columns
            ):
                series = df.set_index("ObsTime")["EUValue"]
                series.name = mnemonic
                combined_dfs.append(series)
            mnemonic_count += 1

        if pbarparent is not None:
            pbarparent.progressbar_increment(
                len(mnemonics), len(mnemonics), "Creating shared time index."
            )

        if combined_dfs:
            combined_df = pd.concat(combined_dfs, axis=1)
            combined_df = combined_df[~combined_df.index.duplicated(keep="first")]
            combined_df.index.name = "ObsTime"
            end_cpu = time.process_time()
            if verbose:
                print('DF combine processing time: ', end_cpu - start_cpu)
            return combined_df
        else:
            return pd.DataFrame()


def save_edb_data(df: pd.DataFrame,
                  output_path: str,
                  file_format: str = 'csv'):
    if file_format == 'csv':
        df.to_csv(output_path, index=False)
    elif file_format == 'parquet':
        df.to_parquet(output_path, index=False)
    elif file_format == 'hdf5':
        df.to_hdf(output_path, key='data', mode='w')
    elif file_format == 'pickle':
        df.to_pickle(output_path)
    else:
        raise ValueError(f"Unsupported format: {file_format}")

    print(f"Saved {len(df)} rows to {output_path}")


# ---------------------------------------------------------------------------
# Group / mnemonic helpers
# ---------------------------------------------------------------------------

def _load_yaml(path: str) -> dict:
    """Load a YAML file, raising a clear error if PyYAML isn't installed."""
    try:
        import yaml  # type: ignore
    except ImportError as e:
        raise ImportError(
            "PyYAML is required to use --plot-groups. "
            "Install it with `pip install pyyaml`."
        ) from e
    with open(path, "r") as f:
        return yaml.safe_load(f)

# --- Add near the other top-level helpers, after _load_yaml() ------------

def _load_derived_definitions(groups_config: str) -> dict:
    """Return the `derived:` mapping from the groups YAML (may be empty)."""
    cfg = _load_yaml(groups_config)
    return cfg.get("derived", {}) or {}


def _resolve_derived_dependencies(
    derived_defs: dict, requested: Sequence[str]
) -> List[str]:
    """
    Given a list of requested (possibly derived) mnemonics, walk the
    `depends_on` graph and return the full ordered list of derived names
    that need to be evaluated (topologically sorted, dependencies first).
    """
    order: List[str] = []
    visiting: set = set()
    visited: set = set()

    def visit(name: str):
        if name in visited or name not in derived_defs:
            return
        if name in visiting:
            raise ValueError(f"Cyclic dependency detected involving '{name}'")
        visiting.add(name)
        for dep in (derived_defs[name].get("depends_on") or []):
            visit(dep)
        visiting.remove(name)
        visited.add(name)
        order.append(name)

    for r in requested:
        visit(r)
    return order


def _collect_raw_inputs(derived_defs: dict, derived_names: Sequence[str]) -> List[str]:
    """
    Return the de-duplicated list of raw MAST mnemonics (WFI_MCU_B_ANALOG_N etc.)
    needed to evaluate the given derived mnemonics.  Any input that is itself
    another derived mnemonic is excluded (it will be computed, not queried).
    """
    raw = []
    seen: set = set()
    for name in derived_names:
        spec = derived_defs.get(name, {}) or {}
        inputs = (spec.get("inputs") or {})
        for _, src in inputs.items():
            if src in derived_defs:
                continue  # dependency handled via depends_on
            if src not in seen:
                seen.add(src)
                raw.append(src)
    return raw


def _partition_derived(
    mnemonics: Sequence[str], derived_defs: dict
) -> Tuple[List[str], List[str]]:
    """
    Split a mnemonic list into (raw_mnemonics, derived_mnemonics)
    based on the `derived:` section of the YAML.
    """
    derived, raw = [], []
    for m in mnemonics:
        (derived if m in derived_defs else raw).append(m)
    return raw, derived


def _compute_derived_columns(
    combined_df: "pd.DataFrame",
    derived_defs: dict,
    derived_names: Sequence[str],
    verbose: bool = True,
) -> "pd.DataFrame":
    """
    Evaluate derived mnemonic formulas and add them as columns to a wide
    (combined) DataFrame indexed by ObsTime.  Derived mnemonics are computed
    in dependency order (`depends_on`), so intermediates like
    WFI_MCU_B_BOARD_T_RAW become available for downstream formulas.
    """
    import numpy as np

    if combined_df is None or combined_df.empty:
        return combined_df

    # Topologically ordered list
    order = _resolve_derived_dependencies(derived_defs, derived_names)

    for name in order:
        spec = derived_defs[name]
        inputs = (spec.get("inputs") or {})
        formula = spec.get("formula")
        if not formula:
            if verbose:
                print(f"⚠ Derived '{name}' has no formula; skipping.")
            continue

        # Build local namespace: each variable -> pd.Series (or scalar).
        local_ns = {"np": np}
        missing = []
        for var, src in inputs.items():
            if src in combined_df.columns:
                local_ns[var] = combined_df[src]
            else:
                missing.append(src)

        if missing:
            if verbose:
                print(f"⚠ Cannot compute '{name}': missing input columns {missing}")
            combined_df[name] = np.nan
            continue

        try:
            combined_df[name] = eval(  # noqa: S307 (trusted YAML input)
                formula, {"__builtins__": {}}, local_ns
            )
            if verbose:
                units = spec.get("units", "")
                print(f"  ✓ Computed derived mnemonic '{name}' [{units}]")
        except Exception as e:
            if verbose:
                print(f"✗ Failed to evaluate '{name}': {e}")
            combined_df[name] = np.nan

    return combined_df
def list_available_groups(groups_config: str) -> List[Tuple[str, str, int]]:
    """
    Return a list of (group_name, label, n_mnemonics) tuples from a YAML config.

    Expected YAML structure:
        groups:
          <group_name>:
            label: <str>
            mnemonics: [<mnemonic>, ...]
            y_label: <str>   # optional
    """
    cfg = _load_yaml(groups_config)
    groups = cfg.get("groups", {}) or {}
    out = []
    for name, info in groups.items():
        label = (info or {}).get("label", name)
        mnems = (info or {}).get("mnemonics", []) or []
        out.append((name, label, len(mnems)))
    return out


def _mnemonics_of(group_info: dict) -> List[str]:
    """Extract mnemonic *names* from a group entry, accepting list or dict form."""
    m = (group_info or {}).get("mnemonics") or {}
    if isinstance(m, dict):
        return list(m.keys())
    if isinstance(m, list):
        return list(m)
    return []


def _labels_of(group_info: dict) -> dict:
    """Extract {mnemonic: description} from a group entry (empty if list form)."""
    m = (group_info or {}).get("mnemonics") or {}
    if isinstance(m, dict):
        return dict(m)
    return {}


def extract_mnemonics_from_groups(
    groups_config: str,
    selected_groups: Optional[Sequence[str]] = None,
    verbose: bool = True,
) -> List[str]:
    cfg = _load_yaml(groups_config)
    groups = cfg.get("groups", {}) or {}
    selected_set = set(selected_groups) if selected_groups else None

    mnemonics: List[str] = []
    seen = set()
    for name, info in groups.items():
        if selected_set is not None and name not in selected_set:
            continue
        for mn in _mnemonics_of(info):
            if mn not in seen:
                seen.add(mn)
                mnemonics.append(mn)
    return mnemonics


def extract_labels_from_groups(
    groups_config: str,
    selected_groups: Optional[Sequence[str]] = None,
) -> dict:
    """
    Return {mnemonic: 'human description'} for all mnemonics in the
    given (or all) groups.  Missing descriptions are omitted.
    """
    cfg = _load_yaml(groups_config)
    groups = cfg.get("groups", {}) or {}
    selected_set = set(selected_groups) if selected_groups else None

    labels: dict = {}
    for name, info in groups.items():
        if selected_set is not None and name not in selected_set:
            continue
        for mn, desc in _labels_of(info).items():
            if desc and mn not in labels:
                labels[mn] = str(desc)
    return labels

def expand_to_raw_mnemonics(
    mnemonics: Sequence[str],
    groups_config: Optional[str],
    verbose: bool = True,
) -> Tuple[List[str], List[str]]:
    """
    Split `mnemonics` into (raw_to_query, derived_to_compute).  Also expands
    each derived mnemonic into its raw ANALOG input channels and includes any
    dependency chain (e.g. WFI_MCU_B_BOARD_T_RAW).

    Returns
    -------
    raw_query_list : list[str]
        The de-duplicated set of mnemonics to actually request from MAST.
    derived_list : list[str]
        The derived mnemonics that should be computed client-side, in
        dependency order.
    """
    if not groups_config:
        return list(mnemonics), []

    derived_defs = _load_derived_definitions(groups_config)
    if not derived_defs:
        return list(mnemonics), []

    raw, derived = _partition_derived(mnemonics, derived_defs)

    # Walk dependencies to include intermediates (e.g. BOARD_T_RAW)
    ordered_derived = _resolve_derived_dependencies(derived_defs, derived)

    # Gather all raw MAST inputs needed
    needed_raw = _collect_raw_inputs(derived_defs, ordered_derived)

    # Merge with any raw mnemonics the user already requested
    seen = set()
    raw_query_list: List[str] = []
    for m in raw + needed_raw:
        if m not in seen:
            seen.add(m)
            raw_query_list.append(m)

    if verbose and ordered_derived:
        print(
            f"  ↳ {len(ordered_derived)} derived mnemonic(s) will be computed "
            f"from {len(needed_raw)} raw MCU-B analog channel(s)."
        )

    return raw_query_list, ordered_derived

def _read_mnemonics_file(path: str) -> List[str]:
    """Read a whitespace/newline/comma-separated list of mnemonics from a file."""
    with open(path, "r") as f:
        text = f.read()
    tokens = []
    for line in text.splitlines():
        # Strip comments starting with '#'
        line = line.split("#", 1)[0]
        for tok in line.replace(",", " ").split():
            tok = tok.strip()
            if tok:
                tokens.append(tok)
    # De-dup, preserve order
    seen, out = set(), []
    for t in tokens:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


# ---------------------------------------------------------------------------
# Python API
# ---------------------------------------------------------------------------

def query_telemetry(
    mnemonics=None, start_time=None, end_time=None, *,
    groups_config=None, selected_groups=None,
    server="mast", api_token=None,
    combine=True, max_retries=3, retry_delay=1.0, verbose=True,
):
    if start_time is None or end_time is None:
        raise ValueError("start_time and end_time are required")

    # Resolve mnemonic list from groups + explicit args
    resolved: List[str] = list(mnemonics) if mnemonics else []
    if groups_config:
        extracted = extract_mnemonics_from_groups(
            groups_config, selected_groups=selected_groups, verbose=verbose
        )
        seen = set(extracted)
        merged = list(extracted)
        for m in resolved:
            if m not in seen:
                seen.add(m)
                merged.append(m)
        resolved = merged

    if not resolved:
        raise ValueError("No mnemonics to query.")

    # NEW: split raw / derived
    raw_to_query, derived_to_compute = expand_to_raw_mnemonics(
        resolved, groups_config, verbose=verbose
    )

    querier = MASTEngDBQuery(api_token=api_token, server=server)
    raw_results = querier.query_multiple_mnemonics(
        raw_to_query, start_time, end_time,
        combine=False, max_retries=max_retries,
        retry_delay=retry_delay, verbose=verbose,
    )

    # Build wide DF and compute derived columns
    combined_dfs = []
    for mnem, df in raw_results.items():
        if df is None or df.empty:
            continue
        if "ObsTime" in df.columns and "EUValue" in df.columns:
            s = df.set_index("ObsTime")["EUValue"]
            s.name = mnem
            combined_dfs.append(s)

    if combined_dfs:
        wide = pd.concat(combined_dfs, axis=1)
        wide = wide[~wide.index.duplicated(keep="first")]
        wide.index.name = "ObsTime"
    else:
        wide = pd.DataFrame()

    if derived_to_compute and not wide.empty and groups_config:
        derived_defs = _load_derived_definitions(groups_config)
        wide = _compute_derived_columns(
            wide, derived_defs, derived_to_compute, verbose=verbose
        )

    if combine:
        return wide

    # combine=False: return dict.  Add synthetic long-form DataFrames for
    # derived mnemonics so callers see them alongside raw results.
    if derived_to_compute and not wide.empty:
        for name in derived_to_compute:
            if name in wide.columns:
                s = wide[name].dropna()
                if not s.empty:
                    raw_results[name] = pd.DataFrame({
                        "ObsTime": s.index,
                        "EUValue": s.values,
                        "mnemonic": name,
                    })
    return raw_results

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

_EPILOG = """\
Examples
--------
# Quickest form — group name(s) from tlm_groups.yaml, last 24 h, interactive plot
roman-telem rcs_pd
roman-telem bank1_leds bank2_leds
roman-telem rcs_pd --days 7

# Save to CSV instead of plotting
roman-telem rcs_pd --output rcs_pd.csv
roman-telem rcs_pd --days 3 --output rcs_pd_3d.csv

# Query specific mnemonics, explicit time range
roman-telem -m WFI_MCE_SRCS_PD1_V WFI_MCE_SRCS_PD2_V \\
            -s "2026-05-01" -e "2026-07-02" --server int \\
            --output data.csv

# One subplot per mnemonic, interactive
roman-telem -m WFI_MCE_SRCS_PD1_V WFI_MCE_SRCS_PD2_V \\
            -s "2026-09-11" -e "2026-09-14" \\
            --plot-per-mnemonic --show

# Discover available groups
roman-telem --list-groups
roman-telem --plot-groups wfi_mce_srcs.yaml --list-groups

# Use a custom groups file
roman-telem --plot-groups wfi_mce_srcs.yaml --select-groups pd1,pd2 \\
            -s "2026-05-01" -e "2026-07-02" --plot-output pd_trending.png
"""


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="roman-telem",
        description=(
            "Query Roman telemetry from MAST EDB. "
            "Quickest usage: roman-telem <group> [group2 ...] "
            "(looks up groups in tlm_groups.yaml, plots last 24 h interactively). "
            "Use -m for individual mnemonics, --output to save data, "
            "--days to change the look-back window."
        ),
        epilog=_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Positional group names (most convenient shorthand)
    p.add_argument(
        "groups",
        nargs="*",
        default=None,
        metavar="GROUP",
        help="One or more group names from the default groups file "
             "(tlm_groups.yaml in cwd, or $ROMAN_TELEM_GROUPS). "
             "Equivalent to --plot-groups <file> --select-groups <group> --show.",
    )

    # Shorthand: --plot <group> [group2 ...]
    p.add_argument(
        "--plot",
        nargs="+",
        dest="plot_shorthand",
        default=None,
        metavar="GROUP",
        help="Alias for positional GROUP arguments.",
    )

    # Mnemonic sources (any combination is allowed)
    p.add_argument(
        "-m",
        nargs="+",
        dest="mnemonics",
        default=None,
        help="One or more mnemonic names (space-separated).",
    )
    p.add_argument(
        "--mnemonics-file",
        default=None,
        help="Path to a text file with mnemonics (whitespace/comma separated, "
             "'#' comments allowed).",
    )
    p.add_argument(
        "--plot-groups",
        default=None,
        help="Path to a YAML file describing plot groups; mnemonics will be "
             "auto-extracted from it.",
    )
    p.add_argument(
        "--select-groups",
        default=None,
        help="Comma-separated list of group names to include from --plot-groups "
             "(default: all groups).",
    )
    p.add_argument(
        "--list-groups",
        action="store_true",
        help="List available groups in --plot-groups (or the default groups file) "
             "and exit.",
    )

    # Time range
    p.add_argument("-s", "--start", default=None,
                   help="Start time (e.g. '2026-05-01' or ISO datetime). "
                        "If omitted, derived from --days before now.")
    p.add_argument("-e", "--end", default=None,
                   help="End time (e.g. '2026-07-02' or ISO datetime). "
                        "If omitted, defaults to now.")
    p.add_argument("--days", type=float, default=1.0,
                   help="Number of past days to retrieve when -s/--start is not "
                        "specified (default: 1).")

    # Server / auth
    p.add_argument("--server", choices=["mast", "int"], default="mast",
                   help="Which MAST EDB server to query (default: mast).")
    p.add_argument("--api-token", default=None,
                   help="MAST API token. Falls back to MAST_API_TOKEN or "
                        "MAST_API_TOKEN_INT env vars.")

# Query behavior
    p.add_argument("--no-combine", action="store_true",
                   help="Do not concatenate results into a single DataFrame; "
                        "instead save one file per mnemonic (suffix will be "
                        "appended to --output).")
    p.add_argument("--max-retries", type=int, default=3,
                   help="Max retry attempts per mnemonic (default: 3).")
    p.add_argument("--retry-delay", type=float, default=1.0,
                   help="Base retry delay in seconds (linear backoff, default: 1.0).")
    p.add_argument("--quiet", "-q", action="store_true",
                   help="Suppress progress output.")

    # Output
    p.add_argument("--output", "-o", default=None,
                   help="Output data file (extension chooses format: "
                        ".csv, .parquet, .h5/.hdf5, .pkl/.pickle).")
    p.add_argument("--output-format", default=None,
                   choices=["csv", "parquet", "hdf5", "pickle"],
                   help="Override output format (otherwise inferred from --output).")

    # Plotting
    p.add_argument("--plot-output", default=None,
                   help="If set, produce a plot at this path using roman_telem_plot.")
    p.add_argument("--plot-layout", default="vertical",
                   choices=["vertical", "horizontal", "grid"],
                   help="Plot layout when using --plot-groups (default: vertical).")
    p.add_argument("--plot-per-mnemonic", action="store_true",
                   help="Produce one subplot per mnemonic (ignores --plot-groups grouping). "
                        "Requires --plot-output or --show.")
    p.add_argument("--show", action="store_true",
                   help="Display the plot live in a matplotlib window.")

    return p


def _infer_output_format(path: str, override: Optional[str]) -> str:
    if override:
        return override
    lower = path.lower()
    if lower.endswith(".csv"):
        return "csv"
    if lower.endswith(".parquet"):
        return "parquet"
    if lower.endswith(".h5") or lower.endswith(".hdf5"):
        return "hdf5"
    if lower.endswith(".pkl") or lower.endswith(".pickle"):
        return "pickle"
    # Default
    return "csv"


def _handle_list_groups(groups_config: str) -> int:
    """Print available groups from a YAML file and return exit code."""
    try:
        info = list_available_groups(groups_config)
    except FileNotFoundError:
        print(f"ERROR: groups file not found: {groups_config}", file=sys.stderr)
        return 2
    except Exception as e:
        print(f"ERROR: failed to parse {groups_config}: {e}", file=sys.stderr)
        return 2

    if not info:
        print(f"(No groups found in {groups_config})")
        return 0

    # Nicely aligned listing
    name_w = max(len(n) for n, _, _ in info)
    label_w = max(len(l) for _, l, _ in info)
    for name, label, count in info:
        print(f"{name.ljust(name_w)}  -  {label.ljust(label_w)}  "
              f"({count} mnemonic{'s' if count != 1 else ''})")
    print(f"\nTotal: {len(info)} group{'s' if len(info) != 1 else ''}")
    return 0


def _collect_mnemonics_from_args(args: argparse.Namespace,
                                 verbose: bool) -> List[str]:
    """Merge mnemonics from --mnemonics, --mnemonics-file, and --plot-groups."""
    mnemonics: List[str] = []
    seen = set()

    def _add(m: str) -> None:
        if m and m not in seen:
            seen.add(m)
            mnemonics.append(m)

    if args.mnemonics:
        for m in args.mnemonics:
            _add(m)

    if args.mnemonics_file:
        for m in _read_mnemonics_file(args.mnemonics_file):
            _add(m)

    if args.plot_groups:
        selected = None
        if args.select_groups:
            selected = [s.strip() for s in args.select_groups.split(",") if s.strip()]
        for m in extract_mnemonics_from_groups(
            args.plot_groups, selected_groups=selected, verbose=verbose
        ):
            _add(m)

    return mnemonics


def _resolve_default_groups_file() -> Optional[str]:
    """Return the default groups YAML path ($ROMAN_TELEM_GROUPS or tlm_groups.yaml in cwd)."""
    env = os.environ.get("ROMAN_TELEM_GROUPS")
    if env and os.path.isfile(env):
        return env
    local = os.path.join(os.getcwd(), "tlm_groups.yaml")
    if os.path.isfile(local):
        return local
    return None


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    verbose = not args.quiet

    # --- Resolve default groups file when --select-groups used without --plot-groups ---
    if args.select_groups and not args.plot_groups and not args.plot_shorthand:
        args.plot_groups = _resolve_default_groups_file()
        if not args.plot_groups:
            parser.error(
                "--select-groups requires --plot-groups or a default tlm_groups.yaml "
                "in the current directory (or $ROMAN_TELEM_GROUPS)."
            )

    # --- Expand positional groups / --plot shorthand ---
    _group_shorthand = args.plot_shorthand or args.groups or None
    if _group_shorthand:
        if not args.plot_groups:
            args.plot_groups = _resolve_default_groups_file()
        if not args.plot_groups:
            parser.error(
                "Group shorthand requires a groups file. Place tlm_groups.yaml in the "
                "current directory or set $ROMAN_TELEM_GROUPS."
            )
        args.select_groups = ",".join(_group_shorthand)
        args.show = True

    # --- Handle --list-groups (does not require start/end) ---
    if args.list_groups:
        groups_file = args.plot_groups or _resolve_default_groups_file()
        if not groups_file:
            parser.error(
                "--list-groups requires --plot-groups or a default tlm_groups.yaml "
                "in the current directory (or $ROMAN_TELEM_GROUPS)."
            )
        return _handle_list_groups(groups_file)

    # --- Resolve time range ---
    from datetime import timezone
    now = datetime.now(timezone.utc).replace(tzinfo=None)
    if not args.end:
        args.end = now.strftime("%Y-%m-%dT%H:%M:%S")
    if not args.start:
        args.start = (now - pd.Timedelta(days=args.days)).strftime("%Y-%m-%dT%H:%M:%S")
        if verbose:
            print(f"No start time specified — retrieving last {args.days:g} day(s): "
                  f"{args.start} to {args.end} UTC")

    # --- Resolve mnemonic list ---
    try:
        mnemonics = _collect_mnemonics_from_args(args, verbose=verbose)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    except Exception as e:
        print(f"ERROR resolving mnemonics: {e}", file=sys.stderr)
        return 2

    if not mnemonics:
        parser.error("No mnemonics specified. Use --mnemonics, --mnemonics-file, "
                     "or --plot-groups.")

    # NEW: separate raw-vs-derived using the groups YAML `derived:` section
    groups_file = args.plot_groups or _resolve_default_groups_file()
    raw_to_query, derived_to_compute = expand_to_raw_mnemonics(
        mnemonics, groups_file, verbose=verbose
    )

    if verbose:
        print(
            f"Resolved {len(raw_to_query)} raw mnemonic(s) to query "
            f"(+{len(derived_to_compute)} derived to compute)."
        )

    # --- Query MAST EDB ---
    try:
        querier = MASTEngDBQuery(api_token=args.api_token, server=args.server)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    combine = not args.no_combine
    try:
        raw_results = querier.query_multiple_mnemonics(
            raw_to_query,        # <-- pass raw list, not `mnemonics`
            args.start,
            args.end,
            combine=False,
            max_retries=args.max_retries,
            retry_delay=args.retry_delay,
            verbose=verbose,
        )
    except Exception as e:
        print(f"ERROR during query: {e}", file=sys.stderr)
        return 1

    # Build the wide combined DataFrame first (needed to compute derived cols).
    plot_df = pd.concat(
        [df for df in raw_results.values() if df is not None and not df.empty],
        ignore_index=True,
    ) if raw_results else pd.DataFrame()

    combined_dfs = []
    for mnem, df in raw_results.items():
        if df is None or df.empty:
            continue
        if "ObsTime" in df.columns and "EUValue" in df.columns:
            series = df.set_index("ObsTime")["EUValue"]
            series.name = mnem
            combined_dfs.append(series)

    if combined_dfs:
        wide = pd.concat(combined_dfs, axis=1)
        wide = wide[~wide.index.duplicated(keep="first")]
        wide.index.name = "ObsTime"
    else:
        wide = pd.DataFrame()

    # NEW: compute derived columns from the wide DataFrame
    if derived_to_compute and not wide.empty:
        derived_defs = _load_derived_definitions(groups_file)
        wide = _compute_derived_columns(
            wide, derived_defs, derived_to_compute, verbose=verbose
        )

        # Also add derived data back into the long-form plot_df so the
        # existing plotting code (which expects ObsTime/EUValue/mnemonic)
        # can visualize derived quantities too.
        long_rows = []
        for name in derived_to_compute:
            if name in wide.columns:
                s = wide[name].dropna()
                if not s.empty:
                    long_rows.append(
                        pd.DataFrame({
                            "ObsTime": s.index,
                            "EUValue": s.values,
                            "mnemonic": name,
                        })
                    )
        if long_rows:
            plot_df = pd.concat([plot_df, *long_rows], ignore_index=True)

    result = wide if combine else raw_results

    # --- Build label map for column renaming ---
    _label_map: dict = {}
    if args.plot_groups:
        selected_for_labels = (
            [s.strip() for s in args.select_groups.split(",") if s.strip()]
            if args.select_groups else None
        )
        try:
            _label_map = extract_labels_from_groups(
                args.plot_groups, selected_groups=selected_for_labels
            )
        except Exception:
            pass

    # --- Save output ---
    if args.output:
        fmt = _infer_output_format(args.output, args.output_format)
        if combine:
            if isinstance(result, dict):
                result = pd.concat(result.values(), ignore_index=True)
            save_df = result
            if _label_map and hasattr(save_df, "columns"):
                save_df = save_df.rename(columns={
                    m: f"{m} - {desc}" for m, desc in _label_map.items()
                    if m in save_df.columns
                })
            save_edb_data(save_df, args.output, file_format=fmt)
        else:
            # One file per mnemonic
            base, ext = os.path.splitext(args.output)
            if not ext:
                ext = "." + ("csv" if fmt == "csv" else fmt)
            for mnem, df in (result or {}).items():
                if df is None or df.empty:
                    if verbose:
                        print(f"  (skipping empty result for {mnem})")
                    continue
                out_path = f"{base}_{mnem}{ext}"
                save_edb_data(df, out_path, file_format=fmt)

    # --- Display to terminal if no output file ---
    if not args.output:
        if isinstance(result, dict):
            for mnem, df in result.items():
                if df is not None and not df.empty:
                    print(f"\n{'='*60}\n{mnem}\n{'='*60}")
                    print(df)
        else:
            print(result)

    # --- Optional plotting ---
    if args.plot_output or args.show:
        try:
            from roman_telem_plot import plot_telemetry, plot_per_mnemonic  # lazy import
        except ImportError as e:
            print(f"ERROR: unable to import roman_telem_plot: {e}", file=sys.stderr)
            return 1

        selected = None
        if args.select_groups:
            selected = [s.strip() for s in args.select_groups.split(",") if s.strip()]

        label_map = {}
        if args.plot_groups:
            try:
                label_map = extract_labels_from_groups(
                    args.plot_groups, selected_groups=selected
                )
            except Exception as e:
                if verbose:
                    print(f"⚠ Could not load legend labels: {e}")

        if plot_df.empty:
            print("WARNING: no data returned — skipping plot.", file=sys.stderr)
        elif args.plot_per_mnemonic:
            plot_per_mnemonic(
                plot_df,
                mnemonics=mnemonics,
                layout=args.plot_layout,
                output=args.plot_output,
                show=args.show,
                label_map=label_map,
            )
        else:
            plot_telemetry(
                plot_df,
                groups_config=args.plot_groups,
                selected_groups=selected,
                layout=args.plot_layout,
                output=args.plot_output,
                show=args.show,
                label_map=label_map,
            )

        if verbose and args.plot_output:
            print(f"Wrote plot to {args.plot_output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())