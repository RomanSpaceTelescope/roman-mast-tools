"""
roman_telem_mast.py
-------------------
Bridge between roman_telem (time-series) and roman_mast (exposure metadata).

Given a time range, query MAST for every Roman exposure that overlaps that
window and group them into contiguous chunks per PROGRAM.  Returns span dicts
suitable for shading behind telemetry plots.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional
import pandas as pd

import roman_mast
from roman_mast import list_data, parse_visit_id, Exposure


@dataclass
class ProgramSpan:
    """One contiguous run of exposures within a single APT program."""
    program: int
    start: pd.Timestamp                       # earliest exposure_start_time
    end: pd.Timestamp                         # latest exposure_end_time
    n_exposures: int
    visit_ids: List[str] = field(default_factory=list)

    @property
    def key(self) -> int:
        return self.program

    @property
    def label(self) -> str:
        return f"Program {self.program}"


def _exposure_time_window(exp: Exposure) -> Optional[tuple]:
    """Best-effort (start, end) as pandas Timestamps.  Returns None if unknown."""
    t0 = getattr(exp, "exposure_start_time", None)
    t1 = (getattr(exp, "exposure_end_time", None)
          or getattr(exp, "exposure_stop_time", None)
          or getattr(exp, "exposure_mid_time", None))
    if t0 is None:
        return None
    try:
        t0 = pd.Timestamp(str(t0))
    except Exception:
        return None
    if t1 is not None:
        try:
            t1 = pd.Timestamp(str(t1))
        except Exception:
            t1 = t0
    else:
        # Fall back to a nominal 3-minute exposure if only start is known.
        t1 = t0 + pd.Timedelta(minutes=3)
    return t0, t1


def query_program_spans(
    start_time,
    end_time,
    *,
    gap_minutes: float = 30.0,
    program: Optional[int] = None,
    detector: Optional[str] = "WFI04",   # any one SCA is enough to tag a chunk
    kinds: Optional[str] = "all",        # product kinds to include (all → includes darks)
    token: Optional[str] = None,
    server: Optional[str] = None,
    verbose: bool = True,
) -> List[ProgramSpan]:
    """Query MAST and return contiguous ProgramSpans overlapping [start_time, end_time].

    A "contiguous chunk" is a run of exposures in the same PROGRAM whose
    successive inter-exposure gap is less than gap_minutes.  Execution plan,
    pass, segment, etc. are not considered for chunk boundaries.

    Parameters
    ----------
    kinds : str or None
        Product kinds to query: 'all' (default, includes L1 uncals and darks),
        or a comma-separated string like 'cal,uncal'. Only affects timing
        metadata collection; no files are downloaded. Set to 'all' to ensure
        pure-dark programs (which never produce L2 cals) still appear.

    Notes
    -----
    * By default we query a single detector (WFI04) to keep the row count
      down; every SCA of a given exposure has the same timing metadata, so
      one is sufficient for span construction.
    * The MAST query is narrowed by `program` if given; otherwise all Roman
      exposures in the underlying archive are pulled and then filtered
      locally against [start_time, end_time].
    * Using kinds='all' ensures dark programs (which only produce L1 uncals)
      are included in the span overlay.
    """
    if verbose:
        print(f"[roman_telem_mast] Querying MAST for exposures "
              f"({start_time} .. {end_time}, kinds={kinds})")

    res = list_data(
        program=program,
        detector=detector,
        kinds=kinds,
        data_level=None,  # let `kinds` drive selection, not data_level
        token=token,
        server=server,
    )

    if verbose:
        print(f"[roman_telem_mast] MAST returned {res.n_results} rows, "
              f"{res.n_exposures} exposures across {res.n_products} files.")

    t_lo = pd.Timestamp(str(start_time))
    t_hi = pd.Timestamp(str(end_time))

    rows: List[tuple] = []   # (program, t0, t1, visit_id)
    for exp in res.exposures:
        win = _exposure_time_window(exp)
        if win is None:
            continue
        t0, t1 = win
        # Overlap with requested window
        if t1 < t_lo or t0 > t_hi:
            continue
        try:
            prog = parse_visit_id(exp.visit_id)["program"]
        except Exception:
            continue
        rows.append((prog, t0, t1, exp.visit_id))

    if not rows:
        if verbose:
            print("[roman_telem_mast] No exposures overlap the requested window.")
        return []

    # Sort by (program, t0) so gap detection is straightforward.
    rows.sort(key=lambda r: (r[0], r[1]))

    gap = pd.Timedelta(minutes=gap_minutes)
    spans: List[ProgramSpan] = []
    cur: Optional[ProgramSpan] = None

    for prog, t0, t1, vid in rows:
        if (cur is None
                or cur.program != prog
                or t0 - cur.end > gap):
            if cur is not None:
                spans.append(cur)
            cur = ProgramSpan(program=prog, start=t0, end=t1,
                              n_exposures=1, visit_ids=[vid])
        else:
            cur.end = max(cur.end, t1)
            cur.n_exposures += 1
            if vid not in cur.visit_ids:
                cur.visit_ids.append(vid)
    if cur is not None:
        spans.append(cur)

    # Sort final spans chronologically for cleaner plotting/logging.
    spans.sort(key=lambda s: s.start)

    if verbose:
        print(f"[roman_telem_mast] Built {len(spans)} contiguous span(s) "
              f"from {len(rows)} exposures.")
        for s in spans:
            print(f"  {s.label}  {s.start}  →  {s.end}  "
                  f"({s.n_exposures} exp)")
    return spans
