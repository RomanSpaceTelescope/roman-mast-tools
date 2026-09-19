"""
roman_telem_cars.py
-------------------
Load the CAR (Commissioning Activity Report) summary CSV and expose it as
CARSpan objects with CAR metadata and program assignment.

CSV columns (no header, or with a header — auto-detected):
    CAR_Number, CAR_Name, APT_Program, Start_UTC, Duration

Expected formats:
  - Start_UTC: ISO 8601 format (e.g. "2026-08-30T11:26:00") or blank for unscheduled CARs
  - Duration: "HH:MM:SS", "HHH:MM:SS" (multi-day hour counts), pandas-style "D days HH:MM:SS",
    or Excel serial-date format "1900-01-DD HH:MM:SS" (duration bug)
  - APT_Program: integer or blank (for launch sequence items, etc.)
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import List, Optional

import pandas as pd


@dataclass
class CARSpan:
    """One contiguous CAR activity window, tagged with program."""
    car_number: str                # e.g. "CAR-173.4"
    car_name: str                  # e.g. "Count-rate Dependent Nonlinearity (Direct)"
    program: Optional[int]         # APT program (e.g. 1025) or None if unassigned
    start: pd.Timestamp
    end: pd.Timestamp

    @property
    def key(self) -> Optional[int]:
        return self.program

    @property
    def label(self) -> str:
        if self.program is None:
            return self.car_number
        return f"{self.car_number}/{self.program}"


_HMS_RE = re.compile(r"^\s*(\d+):(\d+):(\d+)\s*$")
_EXCEL_RE = re.compile(r"^\s*1900-01-(\d+)\s+(\d+):(\d+):(\d+)\s*$")


def _parse_duration(s) -> Optional[pd.Timedelta]:
    """Parse a duration string.  Returns None if unparseable/empty.

    Handles:
      - 'HH:MM:SS' (with unbounded H for multi-day durations)
      - '1900-01-DD HH:MM:SS' (Excel's serial-date bug for durations)
      - Anything pd.to_timedelta accepts (e.g. '1 days 02:00:00')

    Excel's serial-date bug: when Excel stores a duration as a date, it uses:
      1900-01-DD HH:MM:SS where DD corresponds to (actual_days + 1) due to a
      leap-year bug. Additionally, there's a 2-hour offset. So:
        1900-01-02 00:00:00 → (2-1) days + 00:00:00 + 2h = 26 hours
    """
    if s is None:
        return None
    if isinstance(s, pd.Timedelta):
        return s
    s = str(s).strip()
    if not s or s.lower() in ("nan", "none"):
        return None

    # Excel serial-date bug: '1900-01-02 00:00:00' really means 26 hours.
    m = _EXCEL_RE.match(s)
    if m:
        dd, hh, mm, ss = map(int, m.groups())
        # Excel duration formula: (DD - 1) days + HH:MM:SS + 2 hour offset
        days = dd - 1
        hours = hh + 2
        return pd.Timedelta(days=days, hours=hours, minutes=mm, seconds=ss)

    # Plain HH:MM:SS (with H possibly > 24)
    m = _HMS_RE.match(s)
    if m:
        h, mm, ss = map(int, m.groups())
        return pd.Timedelta(hours=h, minutes=mm, seconds=ss)

    # Last resort: let pandas try
    try:
        return pd.to_timedelta(s)
    except Exception:
        return None


def _parse_program(v: str) -> Optional[int]:
    """Parse APT program number. Returns None if blank or unparseable."""
    if not v:
        return None
    s = str(v).strip()
    if not s or s.lower() in ("nan", "none", ""):
        return None
    try:
        return int(float(s))
    except ValueError:
        return None


def load_car_spans(
    path: str,
    *,
    start_window: Optional[pd.Timestamp] = None,
    end_window: Optional[pd.Timestamp] = None,
    include_unassigned_program: bool = True,
    tz_naive: bool = True,
    verbose: bool = True,
) -> List[CARSpan]:
    """Load CAR spans from a CSV with UTC start timestamps and durations.

    Rows with blank Start_UTC are silently skipped (planned but not-yet-scheduled CARs).
    Rows with blank APT_Program are kept if `include_unassigned_program` is True
    (program=None; gray color in plots, label = just the CAR number).

    Parameters
    ----------
    path : str
        Path to the CAR summary CSV file.
    start_window, end_window : pd.Timestamp, optional
        If both provided, only return spans that overlap [start_window, end_window].
        Can be tz-naive (treated as UTC) or tz-aware; internally normalized to UTC.
    include_unassigned_program : bool
        If False, skip CARs with blank APT_Program.
    tz_naive : bool
        If True (default), return spans with tz-naive timestamps (matching MAST EDB convention).
        If False, return tz-aware UTC timestamps.
    verbose : bool
        Print status messages.

    Returns
    -------
    List[CARSpan]
        Sorted by start time, with timestamps either tz-naive or tz-aware per tz_naive param.
    """
    if not os.path.exists(path):
        if verbose:
            print(f"[roman_telem_cars] File not found: {path}")
        return []

    # Helper: normalize any timestamp to UTC-aware
    def _to_utc(ts):
        ts = pd.Timestamp(ts)
        return ts.tz_localize("UTC") if ts.tz is None else ts.tz_convert("UTC")

    # Read with no assumed header; detect if first row is header-ish.
    df = pd.read_csv(path, header=None,
                     names=["car_number", "car_name",
                            "apt_program", "start_utc", "duration"])

    # If the first row's start_utc doesn't parse as a timestamp, treat it as a header
    # and drop it.
    try:
        pd.Timestamp(df.iloc[0]["start_utc"])
    except Exception:
        df = df.iloc[1:].reset_index(drop=True)

    spans: List[CARSpan] = []
    n_skipped_no_start = 0
    n_skipped_bad_dur = 0

    for _, row in df.iterrows():
        start_raw = row["start_utc"]
        if not start_raw or (isinstance(start_raw, float) and pd.isna(start_raw)):
            n_skipped_no_start += 1
            continue
        try:
            start = _to_utc(start_raw)
        except Exception:
            n_skipped_no_start += 1
            continue

        dur = _parse_duration(row["duration"])
        if dur is None:
            n_skipped_bad_dur += 1
            continue

        prog = _parse_program(row["apt_program"])
        if prog is None and not include_unassigned_program:
            continue

        spans.append(CARSpan(
            car_number=str(row["car_number"]).strip(),
            car_name=str(row["car_name"]).strip(),
            program=prog,
            start=start,
            end=start + dur,
        ))

    # Optional window filter: normalize filter window to UTC-aware
    if start_window is not None and end_window is not None:
        start_window = _to_utc(start_window)
        end_window = _to_utc(end_window)
        spans = [s for s in spans
                 if not (s.end < start_window or s.start > end_window)]

    spans = sorted(spans, key=lambda s: s.start)

    # Strip timezone if requested (for compatibility with MAST EDB tz-naive convention)
    if tz_naive:
        for s in spans:
            if s.start.tz is not None:
                s.start = s.start.tz_localize(None)
                s.end = s.end.tz_localize(None)

    if verbose:
        print(f"[roman_telem_cars] Loaded {len(spans)} CAR span(s) from {path}")
        if n_skipped_no_start:
            print(f"  ({n_skipped_no_start} rows skipped — no Start_UTC)")
        if n_skipped_bad_dur:
            print(f"  ({n_skipped_bad_dur} rows skipped — bad duration)")
        for s in spans[:5]:
            print(f"  {s.car_number:12} P{s.program} {s.start} → {s.end}  ({s.end - s.start})")
        if len(spans) > 5:
            print(f"  ... and {len(spans) - 5} more")

    return spans
