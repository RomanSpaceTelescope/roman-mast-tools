#!/usr/bin/env python3
"""
roman_telem.py
--------------
Command-line and Python API for retrieving Roman Space Telescope telemetry
(mnemonics) from the MAST Engineering Database, with optional plotting via
`roman_telem_plot`.

Wraps `mast_eng_db_query.MASTEngDBQuery` and adds:
    * Mnemonic input from CLI, text file, or YAML "plot groups" config
    * Group filtering (--select-groups) and discovery (--list-groups)
    * CSV/Parquet/HDF5/Pickle output
    * Optional trending plots via roman_telem_plot.plot_telemetry
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


def extract_mnemonics_from_groups(
    groups_config: str,
    selected_groups: Optional[Sequence[str]] = None,
    verbose: bool = True,
) -> List[str]:
    """
    Extract a de-duplicated list of mnemonics from a YAML plot-groups file.

    Parameters
    ----------
    groups_config : str
        Path to YAML groups file.
    selected_groups : list of str, optional
        If given, only mnemonics from these group names are returned.
        If None, mnemonics from all groups are returned.
    """
    cfg = _load_yaml(groups_config)
    groups = cfg.get("groups", {}) or {}

    selected_set = set(selected_groups) if selected_groups else None
    mnemonics: List[str] = []
    seen = set()

    for name, info in groups.items():
        if selected_set is not None and name not in selected_set:
            if verbose:
                print(f"  ✗ Skipping group: {name}")
            continue
        if verbose:
            print(f"  ✓ Including group: {name}")
        for m in (info or {}).get("mnemonics", []) or []:
            if m not in seen:
                seen.add(m)
                mnemonics.append(m)

    # Warn about any requested group names that were not found
    if selected_set is not None:
        missing = selected_set - set(groups.keys())
        if missing and verbose:
            print(
                f"⚠ Warning: the following selected groups were not found in "
                f"{groups_config}: {sorted(missing)}"
            )

    return mnemonics


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
    mnemonics: Optional[Sequence[str]] = None,
    start_time: Union[str, datetime, None] = None,
    end_time: Union[str, datetime, None] = None,
    *,
    groups_config: Optional[str] = None,
    selected_groups: Optional[Sequence[str]] = None,
    server: str = "mast",
    api_token: Optional[str] = None,
    combine: bool = True,
    max_retries: int = 3,
    retry_delay: float = 1.0,
    verbose: bool = True,
) -> Union[pd.DataFrame, dict]:
    """
    Query Roman telemetry from the MAST Engineering DB.

    You may supply mnemonics directly via ``mnemonics``, or auto-extract them
    from a YAML plot-groups file via ``groups_config`` (optionally filtered by
    ``selected_groups``).

    Returns a combined DataFrame (default) or a dict of DataFrames keyed by
    mnemonic when ``combine=False``.
    """
    if start_time is None or end_time is None:
        raise ValueError("start_time and end_time are required")

    # Resolve mnemonic list
    resolved: List[str] = list(mnemonics) if mnemonics else []
    if groups_config:
        extracted = extract_mnemonics_from_groups(
            groups_config, selected_groups=selected_groups, verbose=verbose
        )
        # Merge (extracted first, then any explicit extras) with de-dup
        seen = set(extracted)
        merged = list(extracted)
        for m in resolved:
            if m not in seen:
                seen.add(m)
                merged.append(m)
        resolved = merged

    if not resolved:
        raise ValueError(
            "No mnemonics to query. Provide `mnemonics` and/or `groups_config`."
        )

    if verbose:
        print(f"Querying {len(resolved)} mnemonic(s) from server='{server}' "
              f"between {start_time} and {end_time}")

    querier = MASTEngDBQuery(api_token=api_token, server=server)
    return querier.query_multiple_mnemonics(
        resolved,
        start_time,
        end_time,
        combine=combine,
        max_retries=max_retries,
        retry_delay=retry_delay,
        verbose=verbose,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

_EPILOG = """\
Examples
--------
# 1. Query a few mnemonics directly
roman-telem --mnemonics WFI_MCE_SRCS_PD1_V WFI_MCE_SRCS_PD2_V \\
            -s "2026-05-01" -e "2026-07-02" --server int \\
            --output data.csv

# 2. Discover available groups in a YAML plot-groups config
roman-telem --plot-groups wfi_itps_groups.yaml --list-groups

# 3. Query a subset of groups (much faster than "all")
roman-telem --plot-groups wfi_itps_groups.yaml \\
            --select-groups "fps_mcu_currents,icdh_voltages" \\
            -s "2026-05-01" -e "2026-07-02" --server int

# 4. Query, save data, and produce a trending plot
roman-telem --plot-groups wfi_itps_groups.yaml \\
            --select-groups "fps_mcu_currents,icdh_voltages" \\
            -s "2026-05-01" -e "2026-07-02" --server int \\
            --output selected_data.csv \\
            --plot-output selected_trending.png
"""


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="roman-telem",
        description="Query (and optionally plot) Roman telemetry from MAST EDB.",
        epilog=_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Mnemonic sources (any combination is allowed)
    p.add_argument(
        "--mnemonics", "-m",
        nargs="+",
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
        help="List available groups in --plot-groups and exit.",
    )

    # Time range
    p.add_argument("-s", "--start", default=None,
                   help="Start time (e.g. '2026-05-01' or ISO datetime).")
    p.add_argument("-e", "--end", default=None,
                   help="End time (e.g. '2026-07-02' or ISO datetime).")

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


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    verbose = not args.quiet

    # --- Handle --list-groups (does not require start/end) ---
    if args.list_groups:
        if not args.plot_groups:
            parser.error("--list-groups requires --plot-groups")
        return _handle_list_groups(args.plot_groups)

    # --- Require start/end for actual queries ---
    if not args.start or not args.end:
        parser.error("--start/-s and --end/-e are required (unless using --list-groups)")

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

    if verbose:
        print(f"Resolved {len(mnemonics)} mnemonic(s) to query.")

    # --- Query MAST EDB ---
    try:
        querier = MASTEngDBQuery(api_token=args.api_token, server=args.server)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    combine = not args.no_combine
    try:
        # When plotting, we need the per-mnemonic long-format DataFrames
        # (each carries a 'mnemonic' column), so query with combine=False
        # and build the combined wide DataFrame ourselves if needed.
        raw_results = querier.query_multiple_mnemonics(
            mnemonics,
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

    plot_df = pd.concat(
        [df for df in raw_results.values() if df is not None and not df.empty],
        ignore_index=True,
    ) if raw_results else pd.DataFrame()

    if combine:
        combined_dfs = []
        for mnem, df in raw_results.items():
            if df is None or df.empty:
                continue
            if "ObsTime" in df.columns and "EUValue" in df.columns:
                series = df.set_index("ObsTime")["EUValue"]
                series.name = mnem
                combined_dfs.append(series)
        if combined_dfs:
            result = pd.concat(combined_dfs, axis=1)
            result = result[~result.index.duplicated(keep="first")]
            result.index.name = "ObsTime"
        else:
            result = pd.DataFrame()
    else:
        result = raw_results

    # --- Save output ---
    if args.output:
        fmt = _infer_output_format(args.output, args.output_format)
        if combine:
            if isinstance(result, dict):
                # Shouldn't happen when combine=True, but guard anyway
                result = pd.concat(result.values(), ignore_index=True)
            save_edb_data(result, args.output, file_format=fmt)
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
            from roman_telem_plot import plot_telemetry  # lazy import
        except ImportError as e:
            print(f"ERROR: unable to import roman_telem_plot: {e}", file=sys.stderr)
            return 1

        selected = None
        if args.select_groups:
            selected = [s.strip() for s in args.select_groups.split(",") if s.strip()]

        try:
            plot_telemetry(
                plot_df,
                groups_config=args.plot_groups,
                selected_groups=selected,
                layout=args.plot_layout,
                output=args.plot_output,
                show=args.show,
            )
        except TypeError:
            # If the installed roman_telem_plot has a different signature, retry
            # with a minimal call so we degrade gracefully.
            plot_telemetry(plot_df, output=args.plot_output, show=args.show)

        if verbose and args.plot_output:
            print(f"Wrote plot to {args.plot_output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())