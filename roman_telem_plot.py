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
    "euvalue", "value", "val", "eng_value", "engineering_value",
    "raw_value", "data", "y",
)


def _apply_smoothing(series: pd.Series, window: int) -> pd.Series:
    """
    Rolling-mean smoothing. window=1 (or <=1, or None) returns the series
    unchanged. Uses min_periods=1 so the ends aren't NaN, and center=True
    so features aren't shifted in time.
    """
    if window is None or window <= 1:
        return series
    return series.rolling(window=int(window), min_periods=1, center=True).mean()


def _group_smoothing(groups_cfg: Dict[str, dict], group_name: str) -> int:
    """Look up the smoothing window for a named group; default 1."""
    try:
        gdef = groups_cfg.get(group_name, {}) or {}
        return int(gdef.get("smoothing", 1))
    except (TypeError, ValueError):
        return 1


def _find_group_for_mnemonic(groups_cfg: Dict[str, dict], mnem: str) -> Optional[str]:
    """Find the group that contains the given mnemonic."""
    for gname, gdef in groups_cfg.items():
        mnems = (gdef or {}).get("mnemonics") or {}
        # mnemonics can be either a list or a dict — handle both
        names = mnems if isinstance(mnems, list) else list(mnems.keys())
        if mnem in names:
            return gname
    return None


def _pretty_label(mnem: str, label_map: Optional[Dict[str, str]]) -> str:
    if not label_map:
        return mnem
    desc = label_map.get(mnem)
    return f"{mnem} - {desc}" if desc else mnem


def _truncate(s: str, n: int = 50) -> str:
    """Truncate a string to at most n characters, with an ellipsis if trimmed."""
    return s if len(s) <= n else s[: n - 1] + "…"


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
    """Load a plot-groups YAML file and return the full config dict."""
    try:
        import yaml  # type: ignore
    except ImportError as e:
        raise ImportError(
            "PyYAML is required to load plot-groups config. "
            "Install with `pip install pyyaml`."
        ) from e
    with open(path, "r") as f:
        cfg = yaml.safe_load(f) or {}
    return cfg


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


def _get_groups_from_config(cfg: dict) -> Dict[str, dict]:
    """Extract the 'groups' section from a full config dict."""
    return cfg.get("groups", {}) or {}


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
    label_map: Optional[Dict[str, str]] = None,
    smoothing: int = 1,
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
        # Map known discrete string states to numeric, then coerce the rest
        _STATE_MAP = {"DIS": 0, "ENA": 1, "OFF": 0, "ON": 1,
                      "FALSE": 0, "TRUE": 1, "LOW": 0, "HIGH": 1}
        yvals = sub[vcol].map(
            lambda v: _STATE_MAP.get(str(v).strip().upper(), v)
        )
        yvals = pd.to_numeric(yvals, errors="coerce")

        if smoothing > 1:
            # Faint raw trace behind the smoothed line
            ax.plot(
                sub[tcol], yvals,
                marker=marker, linestyle=linestyle,
                markersize=markersize, alpha=0.25,
                color=None, label=None,
            )
            yvals_smooth = _apply_smoothing(yvals, smoothing)
            ax.plot(
                sub[tcol], yvals_smooth,
                marker=marker, linestyle=linestyle,
                markersize=markersize, alpha=alpha,
                label=f"{_pretty_label(m, label_map)} (rolling N={smoothing})",
            )
        else:
            ax.plot(
                sub[tcol], yvals,
                marker=marker, linestyle=linestyle,
                markersize=markersize, alpha=alpha,
                label=_pretty_label(m, label_map),
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
    label_map: Optional[Dict[str, str]] = None,
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
        label_map=label_map,
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
    label_map: Optional[Dict[str, str]] = None,
    smoothing_override: Optional[int] = None,
    program_spans: Optional[list] = None,
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
        base_h_per_row = 2.2
        # Add extra height for program-label annotations if present
        # CAR legend will be managed by the dedicated phantom axes
        extra_h = 0.8 if program_spans else 0.2
        figsize = (base_w * (ncols / max(1, min(ncols, 2))),
                   base_h_per_row * nrows + extra_h)

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=figsize, dpi=dpi,
        sharex=sharex, squeeze=False,
        layout="constrained",
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
        smoothing = smoothing_override if smoothing_override is not None \
                    else (int(ginfo.get("smoothing", 1)) if ginfo else 1)

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
            label_map=label_map,
            smoothing=smoothing,
        )
        total_plotted += plotted
        total_missing += missing

    # Hide any leftover unused axes (grid layout with n < nrows*ncols)
    for j in range(n, nrows * ncols):
        axes_flat[j].set_visible(False)

    # Shade program spans on each subplot and add callout annotations
    if program_spans:
        bottom_ax = axes_flat[n - 1] if n > 0 else None
        for ax in axes_flat[:n]:
            _shade_program_spans(ax, program_spans)
        # Rotate dates on bottom axis and hide on others
        if bottom_ax:
            for lbl in bottom_ax.get_xticklabels():
                lbl.set_rotation(45)
                lbl.set_ha("right")
                lbl.set_rotation_mode("anchor")
            # Add DOY to date format
            import matplotlib.dates as mdates
            bottom_ax.xaxis.set_major_formatter(mdates.DateFormatter("%j | %m-%d %H"))
        for ax in axes_flat[:n-1]:
            for lbl in ax.get_xticklabels():
                lbl.set_visible(False)
        # Add callout annotations on the bottom axis
        if bottom_ax:
            bottom_ax.set_xlabel("")
            fig.canvas.draw()
            # Choose label content based on span type
            if hasattr(program_spans[0], "car_name"):
                label_fn = lambda s: _truncate(
                    f"{s.car_number}/{s.program}: {s.car_name}", 50
                )
            else:
                label_fn = lambda s: _truncate(f"P{s.program}", 50)
            _annotate_span_callouts(
                bottom_ax, program_spans,
                label_fn=label_fn,
                fontsize=8,
                max_tiers=20,
            )
    else:
        fig.autofmt_xdate()

    if suptitle:
        fig.suptitle(suptitle, fontsize=13)

    # Constrained layout handles all margins automatically

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
    label_map: Optional[Dict[str, str]] = None,
    smoothing_override: Optional[int] = None,
    program_spans: Optional[list] = None,
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
        cfg = _load_groups_config(groups_config)
        groups = _get_groups_from_config(cfg)
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
            label_map=label_map,
            smoothing_override=smoothing_override,
            program_spans=program_spans,
        )

    # No groups_config: single-axes plot
    return plot_mnemonics(
        df,
        mnemonics=mnemonics,
        time_col=time_col, value_col=value_col,
        output=output, figsize=figsize or (12.0, 5.0), dpi=dpi,
        title=suptitle, show=show,
        label_map=label_map,
    )


# ---------------------------------------------------------------------------
# Small utility
# ---------------------------------------------------------------------------

def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.isdir(parent):
        os.makedirs(parent, exist_ok=True)


# ---------------------------------------------------------------------------
# MAST program-span shading
# ---------------------------------------------------------------------------

def _program_color_map(spans):
    """Return {program_id: rgba} using tab10, stable per program."""
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("tab10")
    out = {}
    for s in spans:
        if s.program not in out:
            out[s.program] = cmap(len(out) % 10)
    return out


def _shade_program_spans(ax, spans, *, alpha=0.15):
    """Overlay axvspan boxes for each ProgramSpan on `ax`.
    Colors are stable per program. No legend entries — program labeling
    is handled by `_annotate_program_axis` on the bottom axis."""
    if not spans:
        return
    color_map = _program_color_map(spans)
    for s in spans:
        ax.axvspan(s.start, s.end, color=color_map[s.program],
                   alpha=alpha, zorder=0)


def _annotate_program_axis(
    ax,
    spans,
    *,
    unique_labels_only: bool = False,
    y_offset: float = -0.25,
    row_spacing: float = 0.09,
    fontsize: float = 8.0,
    padding_px: float = 4.0,
    max_rows: int = 2,
) -> int:
    """Annotate axis with span labels, deconflicting overlaps.

    By default, labels every span. If unique_labels_only=True, only labels
    the first occurrence of each unique label (useful when many spans share
    the same label). All spans are still shaded; only labels are deduplicated.

    Labels centered under their span; when overlapping, labels move to
    additional rows up to max_rows. Returns the number of rows used.

    Call after xlim/data are established and after fig.canvas.draw() for
    accurate text-extent measurement."""
    if not spans:
        return 0

    import matplotlib.dates as mdates
    import matplotlib.transforms as mtransforms

    color_map = _program_color_map(spans)
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)

    xlim_lo, xlim_hi = ax.get_xlim()
    row_right_px = []  # rightmost occupied pixel on each row

    def _to_pixels_x(x_data):
        return ax.transData.transform((x_data, 0))[0]

    # Sort spans by start time for left-to-right greedy placement
    ordered = sorted(spans, key=lambda s: s.start)

    # If unique_labels_only, keep only the first occurrence of each label
    if unique_labels_only:
        seen = set()
        spans_to_label = []
        for s in ordered:
            if s.label in seen:
                continue
            seen.add(s.label)
            spans_to_label.append(s)
    else:
        spans_to_label = ordered

    for s in spans_to_label:
        x0 = max(mdates.date2num(s.start), xlim_lo)
        x1 = min(mdates.date2num(s.end), xlim_hi)
        if x1 <= x0:
            continue
        xc = 0.5 * (x0 + x1)
        xc_px = _to_pixels_x(xc)

        # Draw text off-screen to measure its width
        txt = ax.text(
            xc, y_offset,
            f"P{s.program}",
            transform=trans,
            ha="center", va="top",
            fontsize=fontsize,
            color=color_map[s.program],
            fontweight="bold",
            clip_on=False,
            in_layout=True,
        )
        bbox = txt.get_window_extent(renderer=renderer)
        half_w = 0.5 * bbox.width
        left_px = xc_px - half_w
        right_px = xc_px + half_w

        # Find first row where label doesn't collide
        placed = False
        for r, prev_right in enumerate(row_right_px):
            if left_px >= prev_right + padding_px:
                row_right_px[r] = right_px
                if r > 0:
                    txt.set_y(y_offset - r * row_spacing)
                placed = True
                break

        if not placed:
            if len(row_right_px) >= max_rows:
                txt.set_visible(False)
                continue
            # Start a new row
            row_right_px.append(right_px)
            r = len(row_right_px) - 1
            if r > 0:
                txt.set_y(y_offset - r * row_spacing)

    return len(row_right_px)


def _annotate_span_callouts(
    ax,
    spans,
    *,
    fontsize: float = 8.0,
    rotation: float = 45.0,
    max_tiers: int = 20,
    tier_step_inches: float = 0.20,
    base_offset_inches: float = 0.10,
    stub_lw: float = 0.8,
    padding_px: float = 4.0,
    label_fn=None,
    unique_only: bool = False,
):
    """Draw 45°-angled callouts pointing from each shaded span to a label
    below the axis, using variable stub lengths to deconflict overlaps.
    Auto-grows the tier stack as needed up to max_tiers.

    Returns the total number of tiers used.
    """
    if not spans:
        return 0

    import matplotlib.dates as mdates
    import matplotlib.transforms as mtransforms

    color_map = _program_color_map(spans)

    if label_fn is None:
        label_fn = lambda s: s.label

    # Optionally deduplicate labels (rarely needed for CARs)
    if unique_only:
        seen, filtered = set(), []
        for s in sorted(spans, key=lambda x: x.start):
            k = label_fn(s)
            if k in seen:
                continue
            seen.add(k)
            filtered.append(s)
        ordered = filtered
    else:
        ordered = sorted(spans, key=lambda s: s.start)

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    fig_h_in = fig.get_size_inches()[1]
    ax_pos = ax.get_position()
    ax_h_frac = ax_pos.height

    def _in_to_axfrac(inches):
        return inches / (fig_h_in * ax_h_frac)

    base_y = -_in_to_axfrac(base_offset_inches)
    tier_step = _in_to_axfrac(tier_step_inches)

    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    xlim_lo, xlim_hi = ax.get_xlim()

    # For each tier, track list of (x0, x1) horizontal intervals occupied in pixels
    tier_bands = []

    def _overlaps(new_x0, new_x1, bands, pad):
        """Check if [new_x0, new_x1] overlaps any band in bands (with padding)."""
        for (bx0, bx1) in bands:
            if not (new_x1 + pad <= bx0 or new_x0 - pad >= bx1):
                return True
        return False

    for s in ordered:
        # Center x of the span, clipped to visible axis
        x0 = max(mdates.date2num(s.start), xlim_lo)
        x1 = min(mdates.date2num(s.end), xlim_hi)
        if x1 <= x0:
            continue
        xc = 0.5 * (x0 + x1)

        color = color_map[s.program]
        text = label_fn(s)

        # Provisionally draw at tier 0 to measure rotated text extent
        stub_bottom_y = base_y - 0 * tier_step
        txt = ax.text(
            xc, stub_bottom_y, text,
            transform=trans,
            ha="right", va="bottom",
            rotation=rotation,
            rotation_mode="anchor",
            fontsize=fontsize,
            color=color,
            clip_on=False,
        )
        bbox = txt.get_window_extent(renderer=renderer)
        new_x0, new_x1 = bbox.x0, bbox.x1

        # Find lowest tier where this label doesn't overlap any existing band
        placed_tier = None
        for t, bands in enumerate(tier_bands):
            if not _overlaps(new_x0, new_x1, bands, padding_px):
                placed_tier = t
                break

        # If no tier found, create a new one (up to max_tiers)
        if placed_tier is None:
            if len(tier_bands) >= max_tiers:
                txt.set_visible(False)
                continue
            tier_bands.append([])
            placed_tier = len(tier_bands) - 1

        # Move the text to its final tier
        stub_bottom_y = base_y - placed_tier * tier_step
        txt.set_y(stub_bottom_y)

        # Re-measure at final y and record the band
        bbox = txt.get_window_extent(renderer=renderer)
        tier_bands[placed_tier].append((bbox.x0, bbox.x1))

        # Draw the vertical stub from axis to anchor
        ax.plot(
            [xc, xc],
            [0, stub_bottom_y],
            transform=trans,
            color=color,
            lw=stub_lw,
            alpha=0.6,
            solid_capstyle="butt",
            clip_on=False,
            zorder=5,
        )

        # Tiny tick-mark at the top of the stub
        ax.plot(
            [xc], [0],
            transform=trans,
            marker="v",
            markersize=4,
            color=color,
            clip_on=False,
            zorder=6,
        )

    return len(tier_bands)


def _draw_car_legend(fig, spans, *, fontsize: float = 8.0):
    """Draw a bottom-of-figure legend inside a dedicated axes so
    constrained_layout actually reserves space for it.

    Only renders if the spans have a `car_name` attribute (i.e., are CARSpans).
    """
    if not spans or not hasattr(spans[0], "car_name"):
        return

    color_map = _program_color_map(spans)

    # Collect first CAR name per program, in sorted order
    seen = {}
    for s in spans:
        if s.program not in seen:
            seen[s.program] = s.car_name

    n_lines = len(seen)
    if n_lines == 0:
        return

    # Add a new axes at the bottom of the figure, spanning the width.
    # Height in inches, converted to figure fraction, matched to n_lines.
    line_h_in = 1.6 * fontsize / 72.0        # rough line height in inches
    fig_h_in = fig.get_size_inches()[1]
    axes_h = (n_lines * line_h_in) / fig_h_in

    # Use add_axes with a rect that constrained_layout will re-manage.
    # The trick: make it a real axes but hide the frame/ticks.
    legend_ax = fig.add_axes([0, 0, 1, axes_h])
    legend_ax.set_axis_off()

    for i, (prog, name) in enumerate(sorted(seen.items())):
        # Place each line from top to bottom in axes coords
        y = 1.0 - (i + 0.5) / n_lines
        legend_ax.text(
            0.5, y,
            f"P{prog}: {name}",
            ha="center", va="center",
            fontsize=fontsize,
            color=color_map[prog],
            transform=legend_ax.transAxes,
        )


# ---------------------------------------------------------------------------
# Public API: one subplot per mnemonic
# ---------------------------------------------------------------------------

def plot_per_mnemonic(
    df: pd.DataFrame,
    mnemonics: Optional[Sequence[str]] = None,
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
    label_map: Optional[Dict[str, str]] = None,
    program_spans: Optional[list] = None,
):
    """Produce one subplot per mnemonic in a single figure.

    Parameters
    ----------
    df : pd.DataFrame
        Long-form telemetry DataFrame with a 'mnemonic' column.
    mnemonics : list of str, optional
        Mnemonics to plot. If None, all mnemonics present in `df` are used.
    layout : {'vertical', 'horizontal', 'grid'}
        Subplot arrangement.
    output : str, optional
        If given, save the figure to this path.
    """
    import matplotlib.pyplot as plt

    if df is None or (hasattr(df, "empty") and df.empty):
        raise ValueError("plot_per_mnemonic(): input DataFrame is empty.")

    if "mnemonic" not in df.columns:
        raise ValueError("DataFrame must contain a 'mnemonic' column.")

    if mnemonics is None:
        mnemonics = sorted(df["mnemonic"].dropna().unique().tolist())

    mnemonics = list(mnemonics)
    if not mnemonics:
        raise ValueError("plot_per_mnemonic(): no mnemonics to plot.")

    n = len(mnemonics)
    nrows, ncols = _layout_grid(n, layout)

    if figsize is None:
        base_w = 12.0
        base_h_per_row = 2.2
        # Add extra height for program-label annotations if present
        # CAR legend will be managed by the dedicated phantom axes
        extra_h = 0.8 if program_spans else 0.2
        figsize = (base_w * (ncols / max(1, min(ncols, 2))),
                   base_h_per_row * nrows + extra_h)

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=figsize, dpi=dpi,
        sharex=sharex, squeeze=False,
        layout="constrained",
    )

    axes_flat = [axes[r][c] for r in range(nrows) for c in range(ncols)]

    tcol = time_col or _find_time_column(df)
    vcol = value_col or _find_value_column(df, tcol)

    total_plotted = 0
    total_missing = 0

    for i, m in enumerate(mnemonics):
        ax = axes_flat[i]
        plotted, missing = _plot_mnemonics_on_axes(
            ax, df, [m],
            time_col=tcol, value_col=vcol,
            title=_pretty_label(m, label_map),
            x_label="Time" if (not sharex or i >= (nrows - 1) * ncols) else None,
            show_legend=False,
            label_map=label_map,
        )
        total_plotted += plotted
        total_missing += missing

    for j in range(n, nrows * ncols):
        axes_flat[j].set_visible(False)

    # Shade program spans on each subplot and add callout annotations
    if program_spans:
        bottom_ax = axes_flat[n - 1] if n > 0 else None
        for ax in axes_flat[:n]:
            _shade_program_spans(ax, program_spans)
        # Rotate dates on bottom axis and hide on others
        if bottom_ax:
            for lbl in bottom_ax.get_xticklabels():
                lbl.set_rotation(45)
                lbl.set_ha("right")
                lbl.set_rotation_mode("anchor")
            # Add DOY to date format
            import matplotlib.dates as mdates
            bottom_ax.xaxis.set_major_formatter(mdates.DateFormatter("%j | %m-%d %H"))
        for ax in axes_flat[:n-1]:
            for lbl in ax.get_xticklabels():
                lbl.set_visible(False)
        # Add callout annotations on the bottom axis
        if bottom_ax:
            bottom_ax.set_xlabel("")
            fig.canvas.draw()
            # Choose label content based on span type
            if hasattr(program_spans[0], "car_name"):
                label_fn = lambda s: _truncate(
                    f"{s.car_number}/{s.program}: {s.car_name}", 50
                )
            else:
                label_fn = lambda s: _truncate(f"P{s.program}", 50)
            _annotate_span_callouts(
                bottom_ax, program_spans,
                label_fn=label_fn,
                fontsize=8,
                max_tiers=20,
            )
    else:
        fig.autofmt_xdate()

    if suptitle:
        fig.suptitle(suptitle, fontsize=13)

    # Constrained layout handles all margins automatically

    if output:
        _ensure_parent_dir(output)
        fig.savefig(output, dpi=dpi, bbox_inches="tight")
        print(f"Saved plot to {output} "
              f"({total_plotted} traces, {total_missing} missing)")
    if show:
        plt.show()

    return fig


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
    p.add_argument("--smoothing", type=int, default=None,
                   help="Override rolling-mean window (samples) for all plotted "
                        "groups. 1 disables smoothing.")
    args = p.parse_args(argv)

    # Load data
    df = pd.read_csv(args.input) if args.input.endswith('.csv') else pd.read_pickle(args.input)

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
        smoothing_override=args.smoothing,
    )
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(_main())