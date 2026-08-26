#!/usr/bin/env python3
"""
roman_telem_plot.py
-------------------
Plotting utilities for Roman telemetry DataFrames returned by
`roman_telem.query_telemetry` (or `MASTEngDBQuery.query_multiple_mnemonics`).

Main entry points
-----------------
- plot_telemetry(df, groups_config=None, selected_groups=None, layout='vertical',
                 output=None, ...)
    High-level convenience function used by the `roman-telem` CLI and by
    downstream Python code.

- plot_mnemonics(df, mnemonics=None, ax=None, ...)
    Plot one or more mnemonics on a single Axes (no grouping).

- plot_grouped(df, groups, layout='vertical', ...)
    Plot pre-built groups (dict of group_name -> {label, mnemonics, y_label}).

Expected DataFrame shape
------------------------
A "long" DataFrame with at least:
  * a `mnemonic` column,
  * a time column (auto-detected: 'time', 'datetime', 'date', 'timestamp',
    'obs_time', or the first datetime64 column),
  * a value column (auto-detected: 'value', 'val', 'eng_value', 'raw_value',
    or the first numeric non-time column).
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

# Matplotlib is imported lazily so importing this module doesn't force a
# backend selection.


# ---------------------------------------------------------------------------
# Column detection helpers
# ---------------------------------------------------------------------------

_TIME_COL_CANDIDATES = (
    "time", "datetime", "date", "timestamp", "obs_time", "obstime",
    "sample_time", "t",
)
_VALUE_COL_CANDIDATES = (
    "value", "val", "eng_value", "engineering_value",
    "raw_value", "data", "y",
)


def _find_time_column(df: pd.DataFrame) -> str:
    lower = {c.lower(): c for c in df.columns}
    for cand in _TIME_COL_CANDIDATES:
        if cand in lower:
            return lower[cand]
    # Fallback: any datetime64 dtype
    for c in df.columns:
        if pd.api.types.is_datetime64_any_dtype(df[c]):
            return c
    # Fallback: any column with "time" or "date" in the name
    for c in df.columns:
        cl = c.lower()
        if "time" in cl or "date" in cl:
            return c
    raise ValueError(
        f"Could not find a time column in DataFrame. Columns: {list(df.columns)}"
    )


def _find_value_column(df: pd.DataFrame, time_col: str) -> str:
    lower = {c.lower(): c for c in df.columns}
    for cand in _VALUE_COL_CANDIDATES:
        if cand in lower:
            return lower[cand]
    # Fallback: first numeric column that isn't the time column
    for c in df.columns:
        if c == time_col or c == "mnemonic":
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            return c
    raise ValueError(
        f"Could not find a value column in DataFrame. Columns: {list(df.columns)}"
    )


# ---------------------------------------------------------------------------
# YAML groups loading (compatible with roman_telem)
# ---------------------------------------------------------------------------

def _load_groups_config(path: str) -> Dict[str, dict]:
    """Load a plot-groups YAML file and return the `groups` dict."""
    try:
        import yaml  # type: ignore
    except ImportError as e:
        raise ImportError(
            "PyYAML is required to load plot-groups config. "
            "Install with `pip install pyyaml`."
        ) from e
    with open(path, "r") as f:
        cfg = yaml.safe_load(f) or {}
    return cfg.get("groups", {}) or {}


def _filter_groups(
    groups: Dict[str, dict],
    selected: Optional[Sequence[str]],
) -> Dict[str, dict]:
    if not selected:
        return dict(groups)
    sel = list(selected)
    out = {name: groups[name] for name in sel if name in groups}
    missing = [name for name in sel if name not in groups]
    if missing:
        print(f"⚠ Warning: selected groups not found in config: {missing}")
    return out


# ---------------------------------------------------------------------------
# Core plotting primitives
# ---------------------------------------------------------------------------

def _plot_mnemonics_on_axes(
    ax,
    df: pd.DataFrame,
    mnemonics: Sequence[str],
    *,
    time_col: Optional[str] = None,
    value_col: Optional[str] = None,
    title: Optional[str] = None,
    y_label: Optional[str] = None,
    x_label: Optional[str] = "Time",
    marker: str = ".",
    linestyle: str = "-",
    markersize: float = 3.0,
    alpha: float = 0.9,
    show_legend: bool = True,
    legend_loc: str = "best",
    legend_fontsize: float = 8.0,
) -> Tuple[int, int]:
    """Plot the given mnemonics from a long DataFrame onto a single Axes.

    Returns
    -------
    (n_plotted, n_missing)
    """
    if df is None or df.empty:
        ax.set_title(title or "")
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, color="gray")
        return 0, len(mnemonics)

    if "mnemonic" not in df.columns:
        raise ValueError("DataFrame must contain a 'mnemonic' column.")

    tcol = time_col or _find_time_column(df)
    vcol = value_col or _find_value_column(df, tcol)

    # Ensure time is datetime for nice x-axis
    if not pd.api.types.is_datetime64_any_dtype(df[tcol]):
        try:
            df = df.copy()
            df[tcol] = pd.to_datetime(df[tcol], errors="coerce")
        except Exception:
            pass

    plotted = 0
    missing = 0
    for m in mnemonics:
        sub = df[df["mnemonic"] == m]
        if sub.empty:
            missing += 1
            continue
        # Sort by time so lines don't zig-zag
        sub = sub.sort_values(tcol)
        # Coerce value to numeric where possible
        yvals = pd.to_numeric(sub[vcol], errors="coerce")
        ax.plot(
            sub[tcol], yvals,
            marker=marker, linestyle=linestyle,
            markersize=markersize, alpha=alpha,
            label=m,
        )
        plotted += 1

    if title:
        ax.set_title(title)
    if y_label:
        ax.set_ylabel(y_label)
    if x_label:
        ax.set_xlabel(x_label)
    ax.grid(True, which="both", linestyle=":", alpha=0.4)

    if show_legend and plotted > 0:
        ax.legend(loc=legend_loc, fontsize=legend_fontsize, framealpha=0.85)

    return plotted, missing


# ---------------------------------------------------------------------------
# Public API: single-axes plot
# ---------------------------------------------------------------------------

def plot_mnemonics(
    df: pd.DataFrame,
    mnemonics: Optional[Sequence[str]] = None,
    *,
    ax=None,
    time_col: Optional[str] = None,
    value_col: Optional[str] = None,
    title: Optional[str] = None,
    y_label: Optional[str] = None,
    output: Optional[str] = None,
    figsize: Tuple[float, float] = (12.0, 5.0),
    dpi: int = 120,
    show: bool = False,
):
    """Plot one or more mnemonics on a single Axes.

    If `mnemonics` is None, plot every mnemonic present in `df`.
    If `output` is given, the figure is saved to that path.
    """
    import matplotlib.pyplot as plt

    if mnemonics is None:
        if "mnemonic" not in df.columns:
            raise ValueError("DataFrame must contain a 'mnemonic' column.")
        mnemonics = sorted(df["mnemonic"].dropna().unique().tolist())

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure

    _plot_mnemonics_on_axes(
        ax, df, mnemonics,
        time_col=time_col, value_col=value_col,
        title=title, y_label=y_label,
    )

    fig.autofmt_xdate()
    fig.tight_layout()

    if output:
        _ensure_parent_dir(output)
        fig.savefig(output, dpi=dpi, bbox_inches="tight")
        print(f"Saved plot to {output}")
    if show:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Public API: grouped plot
# ---------------------------------------------------------------------------

def _layout_grid(n: int, layout: str) -> Tuple[int, int]:
    """Compute (nrows, ncols) for the given layout style."""
    if n <= 0:
        return 1, 1
    if layout == "vertical":
        return n, 1
    if layout == "horizontal":
        return 1, n
    if layout == "grid":
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))
        return nrows, ncols
    raise ValueError(f"Unknown layout: {layout!r} "
                     "(expected 'vertical', 'horizontal', or 'grid')")


def plot_grouped(
    df: pd.DataFrame,
    groups: Dict[str, dict],
    *,
    layout: str = "vertical",
    time_col: Optional[str] = None,
    value_col: Optional[str] = None,
    output: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: int = 120,
    suptitle: Optional[str] = None,
    show: bool = False,
    sharex: bool = True,
):
    """Plot each group in its own subplot.

    Parameters
    ----------
    df : pd.DataFrame
        Long-form telemetry DataFrame with a 'mnemonic' column.
    groups : dict
        Mapping of group_name -> {'label': str, 'mnemonics': [..],
                                  'y_label': str (optional)}.
    layout : {'vertical', 'horizontal', 'grid'}
        Subplot arrangement.
    output : str, optional
        If given, save the figure to this path.
    """
    import matplotlib.pyplot as plt

    if not groups:
        raise ValueError("No groups provided to plot_grouped().")

    group_items: List[Tuple[str, dict]] = list(groups.items())
    n = len(group_items)
    nrows, ncols = _layout_grid(n, layout)

    # Reasonable default figure size based on layout
    if figsize is None:
        base_w = 12.0
        base_h_per_row = 3.2
        figsize = (base_w * (ncols / max(1, min(ncols, 2))),
                   base_h_per_row * nrows)

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=figsize, dpi=dpi,
        sharex=sharex, squeeze=False,
    )

    # Flatten axes into a list in row-major order
    axes_flat = [axes[r][c] for r in range(nrows) for c in range(ncols)]

    total_plotted = 0
    total_missing = 0

    for i, (gname, ginfo) in enumerate(group_items):
        ax = axes_flat[i]
        ginfo = ginfo or {}
        label = ginfo.get("label", gname)
        mnems = list(ginfo.get("mnemonics", []) or [])
        y_label = ginfo.get("y_label")

        if not mnems:
            ax.set_title(f"{label} (no mnemonics)")
            ax.text(0.5, 0.5, "No mnemonics defined",
                    ha="center", va="center",
                    transform=ax.transAxes, color="gray")
            continue

        plotted, missing = _plot_mnemonics_on_axes(
            ax, df, mnems,
            time_col=time_col, value_col=value_col,
            title=label, y_label=y_label,
            x_label="Time" if (not sharex or i >= (nrows - 1) * ncols) else None,
        )
        total_plotted += plotted
        total_missing += missing

    # Hide any leftover unused axes (grid layout with n < nrows*ncols)
    for j in range(n, nrows * ncols):
        axes_flat[j].set_visible(False)

    if suptitle:
        fig.suptitle(suptitle, fontsize=13)

    fig.autofmt_xdate()
    fig.tight_layout(rect=(0, 0, 1, 0.97 if suptitle else 1.0))

    if output:
        _ensure_parent_dir(output)
        fig.savefig(output, dpi=dpi, bbox_inches="tight")
        print(f"Saved plot to {output} "
              f"({total_plotted} traces, {total_missing} missing)")
    if show:
        plt.show()

    return fig


# ---------------------------------------------------------------------------
# Public API: high-level convenience (used by CLI and downstream code)
# ---------------------------------------------------------------------------

def plot_telemetry(
    df: pd.DataFrame,
    *,
    groups_config: Optional[str] = None,
    selected_groups: Optional[Sequence[str]] = None,
    layout: str = "vertical",
    output: Optional[str] = None,
    mnemonics: Optional[Sequence[str]] = None,
    time_col: Optional[str] = None,
    value_col: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: int = 120,
    suptitle: Optional[str] = None,
    show: bool = False,
):
    """
    High-level plotting entry point.

    Behavior
    --------
    * If `groups_config` is given, load the YAML and produce a subplot per
      group (optionally filtered by `selected_groups`).
    * Otherwise, plot the provided `mnemonics` (or all mnemonics present in
      `df`) on a single Axes.

    Parameters
    ----------
    df : pd.DataFrame
        Long-form telemetry DataFrame with a 'mnemonic' column and time/value
        columns (auto-detected).
    groups_config : str, optional
        Path to a YAML plot-groups file.
    selected_groups : list of str, optional
        Filter groups from `groups_config` by name.
    layout : {'vertical', 'horizontal', 'grid'}
        Subplot arrangement when using grouped plotting.
    output : str, optional
        If given, save the figure to this path.
    mnemonics : list of str, optional
        When `groups_config` is None, restrict single-axes plot to these.
    figsize, dpi, suptitle, show :
        Standard matplotlib knobs.
    """
    if df is None or (hasattr(df, "empty") and df.empty):
        raise ValueError("plot_telemetry(): input DataFrame is empty.")

    if groups_config:
        groups = _load_groups_config(groups_config)
        groups = _filter_groups(groups, selected_groups)
        if not groups:
            raise ValueError(
                "plot_telemetry(): no groups to plot after filtering "
                f"(selected_groups={list(selected_groups) if selected_groups else None})."
            )
        return plot_grouped(
            df, groups,
            layout=layout,
            time_col=time_col, value_col=value_col,
            output=output, figsize=figsize, dpi=dpi,
            suptitle=suptitle, show=show,
        )

    # No groups_config: single-axes plot
    return plot_mnemonics(
        df,
        mnemonics=mnemonics,
        time_col=time_col, value_col=value_col,
        output=output, figsize=figsize or (12.0, 5.0), dpi=dpi,
        title=suptitle, show=show,
    )


# ---------------------------------------------------------------------------
# Small utility
# ---------------------------------------------------------------------------

def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.isdir(parent):
        os.makedirs(parent, exist_ok=True)


# ---------------------------------------------------------------------------
# Optional standalone CLI: python -m roman_telem_plot data.csv --plot-groups ...
# ---------------------------------------------------------------------------

def _main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    p = argparse.ArgumentParser(
        prog="roman-telem-plot",
        description="Plot Roman telemetry data from a saved file.",
    )
    p.add_argument("input", help="Path to CSV/Parquet/HDF5/Pickle telemetry file.")
    p.add_argument("--plot-groups", default=None,
                   help="YAML plot-groups config; if given, produces one "
                        "subplot per group.")
    p.add_argument("--select-groups", default=None,
                   help="Comma-separated list of group names to include.")
    p.add_argument("--mnemonics", "-m", nargs="+", default=None,
                   help="Restrict single-axes plot to these mnemonics "
                        "(ignored if --plot-groups is used).")
    p.add_argument("--layout", default="vertical",
                   choices=["vertical", "horizontal", "grid"])
    p.add_argument("--output", "-o", default=None,
                   help="Output plot path (e.g. trending.png). If omitted, "
                        "the figure is shown interactively.")
    p.add_argument("--suptitle", default=None)
    p.add_argument("--dpi", type=int, default=120)
    args = p.parse_args(argv)

    # Load data using the same helper as roman_telem
    from mast_eng_db_query import load_edb_data
    df = load_edb_data(args.input)

    selected = None
    if args.select_groups:
        selected = [s.strip() for s in args.select_groups.split(",") if s.strip()]

    plot_telemetry(
        df,
        groups_config=args.plot_groups,
        selected_groups=selected,
        mnemonics=args.mnemonics,
        layout=args.layout,
        output=args.output,
        suptitle=args.suptitle,
        dpi=args.dpi,
        show=(args.output is None),
    )
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(_main())