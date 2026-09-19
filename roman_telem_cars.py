"""
roman_telem_cars.py
-------------------
Load the CAR (Commissioning Activity Report) summary CSV and expose it as
ProgramSpan-compatible objects with added `car_number` / `car_name`.
"""
from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field
from typing import List, Optional

import pandas as pd

# Roman launch timestamp: 2026-08-30T11:26:00 UTC
ROMAN_LAUNCH = pd.Timestamp("2026-08-30T11:26:00", tz="UTC")

# CSV columns: CAR_Number, CAR_Name, APT_Program, MET, Start_Time, Duration
# MET format: "L + DDD:HH:MM:SS" (mission elapsed time since launch)
# Start_Time format: "HH:MM:SS" (wall-clock time on MET-day, often rounded)
# Duration format: "HH:MM:SS" (or empty for instantaneous events)

_MET_RE = re.compile(r"L\s*\+\s*(\d+):(\d+):(\d+):(\d+)")


@dataclass
class CARSpan:
    """One contiguous CAR activity window, tagged with program."""
    car_number: str                # e.g. "CAR-173.4"
    car_name: str                  # e.g. "Count-rate Dependent Nonlinearity (Direct)"
    program: int                   # APT program, e.g. 1025
    start: pd.Timestamp
    end: pd.Timestamp

    @property
    def key(self) -> int:
        return self.program

    @property
    def label(self) -> str:
        return self.car_number


def _parse_hms(s: str) -> Optional[pd.Timedelta]:
    """Parse 'HH:MM:SS' → Timedelta.  Returns None on unparseable."""
    if not s:
        return None
    s = str(s).strip()
    if not s or s.lower() in ("delayed", "skip", "-", ""):
        return None
    # Handle edge case: multiline fields in CSV
    s = s.split("\n")[0].strip()
    if not s:
        return None
    m = re.match(r"^(\d+):(\d+):(\d+)$", s)
    if not m:
        return None
    h, mm, ss = map(int, m.groups())
    return pd.Timedelta(hours=h, minutes=mm, seconds=ss)


def _parse_met(met_str: str, launch_time: pd.Timestamp) -> Optional[pd.Timestamp]:
    """Convert 'L + DDD:HH:MM:SS' relative to a launch timestamp."""
    if not met_str:
        return None
    m = _MET_RE.search(str(met_str))
    if not m:
        return None
    d, h, mm, ss = map(int, m.groups())
    return launch_time + pd.Timedelta(days=d, hours=h, minutes=mm, seconds=ss)


def load_car_spans(
    path: str,
    launch_time: pd.Timestamp = ROMAN_LAUNCH,
    *,
    start_window: Optional[pd.Timestamp] = None,
    end_window: Optional[pd.Timestamp] = None,
    verbose: bool = True,
) -> List[CARSpan]:
    """Load CAR spans from a CSV.  Filter to those overlapping the requested
    time window if `start_window` / `end_window` are given.

    Parameters
    ----------
    path : str
        Path to the CAR summary CSV file.
    launch_time : pd.Timestamp
        Reference time for MET=0. Defaults to ROMAN_LAUNCH.
    start_window, end_window : pd.Timestamp, optional
        If both provided, only return spans that overlap [start_window, end_window].
    verbose : bool
        Print status messages.

    Returns
    -------
    List[CARSpan]
        Sorted by start time.
    """
    if not os.path.exists(path):
        if verbose:
            print(f"[roman_telem_cars] File not found: {path}")
        return []

    spans: List[CARSpan] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Parse MET to get the authoritative start time
            met_t = _parse_met(row.get("MET", ""), launch_time)
            if met_t is None:
                continue

            # Parse duration; if missing or invalid, skip (instantaneous events not useful for shading)
            duration = _parse_hms(row.get("Duration", ""))
            if duration is None or duration == pd.Timedelta(0):
                # Skip instantaneous events or unparseable durations
                continue

            start = met_t
            end = start + duration

            # Parse APT program; skip rows without a program number
            try:
                prog_str = row.get("APT_Program", "").strip()
                if not prog_str:
                    continue
                prog = int(prog_str)
            except (TypeError, ValueError):
                continue

            spans.append(CARSpan(
                car_number=row.get("CAR_Number", "").strip(),
                car_name=row.get("CAR_Name", "").strip(),
                program=prog,
                start=start,
                end=end,
            ))

    # Optional window filter
    if start_window is not None and end_window is not None:
        # Ensure tz-aware comparison
        if start_window.tz is None:
            start_window = start_window.tz_localize("UTC")
        if end_window.tz is None:
            end_window = end_window.tz_localize("UTC")
        spans = [s for s in spans
                 if not (s.end < start_window or s.start > end_window)]

    # Sort by start time
    spans = sorted(spans, key=lambda s: s.start)

    if verbose:
        print(f"[roman_telem_cars] Loaded {len(spans)} CAR span(s) from {path}")
        for s in spans[:5]:
            print(f"  {s.car_number} (P{s.program}): {s.start} → {s.end}")
        if len(spans) > 5:
            print(f"  ... and {len(spans) - 5} more")

    return spans
