"""Export Roman WFI metadata to a CSV spreadsheet.

Builds on the streaming primitives in `roman_mast` (same search + selection
model as `roman_fits.py`): search MAST → pick exposures → stream every SCA
into memory → flatten each datamodel's `meta` tree into dot-notation columns
→ write one CSV row per (exposure, SCA).

Public surface
--------------
    flatten_metadata(obj, prefix='')     → dict[str, scalar]  (dot-notation keys)
    extract_row(af, sca, exposure)       → dict for one SCA (data stats + meta)
    extract_rows(af_dict, exposure, ...) → rows for an already-streamed af_dict
    write_csv(rows, output)              → path to written CSV
    write_metadata_csv(af_dict, exp, ...) → one-shot for a single exposure
    export_csv(res, indices, scas, out)  → path to written CSV (streams + writes)

CLI mirrors roman_fits.py — use the standard --program / --pass / --visit-id
/ etc. filters to find exposures, then --exposures / --scas to pick which
to export. See `python export_metadata_csv.py --help` for examples.
"""

from __future__ import annotations

import csv
import os
import warnings
from typing import Any, Optional

import numpy as np
import roman_datamodels as rdm

from roman_mast import (
    Exposure, DataResults, _log, close_streams, parse_int_spec,
)


# ---------------------------------------------------------------------------
# Metadata flattening
# ---------------------------------------------------------------------------

def flatten_metadata(obj, prefix: str = "", visited=None,
                     max_depth: int = 10, current_depth: int = 0) -> dict:
    """Recursively flatten a Roman datamodel meta tree into dot-notation keys.

    Returns a plain ``{"meta.exposure.type": "WFI_IMAGE", ...}`` mapping.
    Handles DNode / dict / __dict__ objects, skips private attrs, and cycles
    are broken via `id()` tracking. Astropy Time-like objects are stringified.
    """
    if visited is None:
        visited = set()
    if current_depth >= max_depth:
        return {}

    obj_id = id(obj)
    if obj_id in visited:
        return {}
    visited.add(obj_id)

    result: dict = {}

    def _handle_value(key, v):
        if v is None:
            return
        # Scalars — take as-is.
        if isinstance(v, (str, int, float, bool)):
            result[key] = v
            return
        # Lists / tuples: stringify short primitive lists, recurse into short
        # mixed lists, skip large ones.
        if isinstance(v, (list, tuple)):
            if len(v) == 0:
                return
            if all(isinstance(x, (str, int, float, bool, type(None))) for x in v):
                result[key] = str(v)[:200]
            elif len(v) < 20:
                result.update(flatten_metadata(
                    v, key, visited, max_depth, current_depth + 1,
                ))
            return
        # Astropy Time / datetime-like → string.
        if hasattr(v, 'iso') or hasattr(v, 'datetime'):
            result[key] = str(v)
            return
        # Nested dict / DNode / __dict__ carrier — recurse.
        if hasattr(v, '__getitem__') or hasattr(v, '__dict__'):
            result.update(flatten_metadata(
                v, key, visited, max_depth, current_depth + 1,
            ))

    try:
        # dict-like / DNode
        if isinstance(obj, dict) or (
            hasattr(obj, '__getitem__') and hasattr(obj, 'keys')
        ):
            try:
                items = list(obj.items())
            except Exception:
                try:
                    items = [(k, obj[k]) for k in obj.keys()]
                except Exception:
                    return {}
            for k, v in items:
                _handle_value(f"{prefix}.{k}" if prefix else str(k), v)

        elif hasattr(obj, '__dict__'):
            for k, v in obj.__dict__.items():
                if k.startswith('_'):
                    continue
                _handle_value(f"{prefix}.{k}" if prefix else str(k), v)

    except Exception:
        pass

    return result


# ---------------------------------------------------------------------------
# Row extraction
# ---------------------------------------------------------------------------

# Priority column order — anything not listed here falls back to alphabetical.
_PRIORITY_KEYS = [
    'visit_id', 'exposure', 'sca',
    'meta.pointing.target_ra',
    'meta.pointing.target_dec',
    'meta.exposure.type',
    'meta.exposure.start_time',
    'meta.exposure.end_time',
    'meta.exposure.frame_time',
    'meta.exposure.effective_exposure_time',
    'meta.pointing.target_aperture',
    'meta.exposure.ma_table_name',
    'meta.statistics.good_pixel_fraction',
    'meta.statistics.image_median',
    'meta.statistics.image_rms',
    'meta.statistics.zodiacal_light',
    'meta.exposure.hga_move',
    'meta.file_date',
    'meta.rcs.active',
    'meta.rcs.bank',
    'meta.rcs.counts',
    'meta.rcs.electronics',
    'meta.rcs.led',
    'meta.observation.execution_plan',
    'meta.observation.exposure',
    'meta.observation.observation',
    'meta.observation.observation_id',
    'meta.observation.pass',
    'meta.observation.program',
    'meta.observation.segment',
    'meta.observation.visit',
    'meta.observation.visit_file_activity',
    'meta.observation.visit_file_group',
    'meta.observation.visit_file_sequence',
    'meta.visit.nexposures',
    'data_shape', 'data_dtype',
    'data_min', 'data_max', 'data_mean',
    'data_nan_pixels', 'data_valid_pixels',
    'meta.source_catalog.tweakreg_catalog_name',
]


def _open_dm(af):
    """Convert an AsdfFile → roman datamodel. Passthrough if already one."""
    if hasattr(af, 'meta') and hasattr(af, 'data'):
        return af
    return rdm.open(af)


def extract_row(af, sca_num: int, exposure: Exposure) -> dict:
    """Extract a single CSV row (data stats + flattened meta) from one SCA."""
    dm = _open_dm(af)

    row: dict = {
        'visit_id': exposure.visit_id,
        'exposure': exposure.exposure,
        'sca': sca_num,
    }

    # Data stats — cheap enough, and handy in the spreadsheet.
    try:
        data = np.asarray(dm.data[...])
        row.update({
            'data_shape':        str(dm.data.shape),
            'data_dtype':        str(dm.data.dtype),
            'data_min':          float(np.nanmin(data)),
            'data_max':          float(np.nanmax(data)),
            'data_mean':         float(np.nanmean(data)),
            'data_valid_pixels': int(np.isfinite(data).sum()),
            'data_nan_pixels':   int(np.isnan(data).sum()),
        })
    except Exception as e:
        row['data_error'] = str(e)

    if hasattr(dm, 'meta'):
        row.update(flatten_metadata(dm.meta, prefix='meta'))

    return row


# ---------------------------------------------------------------------------
# Row extraction from an already-streamed exposure
# ---------------------------------------------------------------------------

def extract_rows(af_dict: dict, exposure: Exposure) -> list:
    """Extract CSV rows from an in-memory `{sca: AsdfFile}` dict.

    Companion to `roman_fits.to_fits_files` / `to_ds9` — those consume the
    same af_dict, so callers who already streamed the exposure for FITS/DS9
    output can drop the metadata CSV for free without re-streaming.
    """
    rows: list = []
    for sca_num in sorted(af_dict):
        af = af_dict[sca_num]
        if af is None:
            _log(f"  SCA {sca_num:02d}: no data — skipping (metadata)")
            continue
        try:
            rows.append(extract_row(af, sca_num, exposure))
        except Exception as e:
            _log(f"  SCA {sca_num:02d}: extract failed: {e}")
            rows.append({
                'visit_id': exposure.visit_id,
                'exposure': exposure.exposure,
                'sca': sca_num,
                'error': str(e),
            })
    return rows


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------

def _order_keys(all_keys):
    """Return the CSV column order — priority keys first, rest alphabetical."""
    all_keys = set(all_keys)
    ordered = [k for k in _PRIORITY_KEYS if k in all_keys]
    remaining = sorted(all_keys - set(ordered))
    return ordered + remaining


def _default_output_for_exposures(exposures):
    """Auto-generate a CSV filename from a list of Exposure objects."""
    if not exposures:
        return 'metadata.csv'
    if len(exposures) == 1:
        exp = exposures[0]
        return f'metadata_{exp.visit_id}_exp{exp.exposure:02d}.csv'
    first, last = exposures[0], exposures[-1]
    if first.visit_id == last.visit_id:
        return (f'metadata_{first.visit_id}'
                f'_exp{first.exposure:02d}-{last.exposure:02d}.csv')
    return f'metadata_{first.visit_id}_plus{len(exposures) - 1}more.csv'


def _default_output(res: DataResults, exposures):
    """Auto-generate a CSV filename from the exposures being exported."""
    return _default_output_for_exposures(exposures)


def write_csv(rows: list, output: str, *, quiet: bool = False) -> str:
    """Write pre-extracted rows to a CSV using the standard column ordering.

    Column order matches `export_csv`: priority keys first, remainder
    alphabetical. Returns the output path.
    """
    if not rows:
        raise RuntimeError("No metadata rows — nothing to write.")

    all_keys: set = set()
    for r in rows:
        all_keys.update(r.keys())
    ordered_keys = _order_keys(all_keys)

    _log(f"Writing CSV → {output}  "
         f"({len(rows)} row(s), {len(ordered_keys)} column(s))")

    with open(output, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=ordered_keys, restval='')
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    if not quiet:
        print(f"\n✓ Wrote {output}")
        print(f"   {len(rows)} row(s), {len(ordered_keys)} column(s)")
    return output


def write_metadata_csv(
    af_dict: dict,
    exposure: Exposure,
    *,
    output: Optional[str] = None,
    quiet: bool = False,
) -> str:
    """One-shot: extract rows from an already-streamed exposure and write CSV.

    Companion to `roman_fits.to_fits_files` / `to_ds9`. Both take the same
    `af_dict`, so a single stream can feed all three outputs.
    """
    if output is None:
        output = _default_output_for_exposures([exposure])
    rows = extract_rows(af_dict, exposure)
    return write_csv(rows, output, quiet=quiet)


def export_csv(
    res: DataResults,
    indices,
    *,
    scas=None,
    output: Optional[str] = None,
    show_progress: bool = True,
) -> str:
    """Stream the selected exposures and write one CSV row per (exp, SCA).

    Parameters
    ----------
    res : DataResults
        Result of `roman_mast.list_data(...)`.
    indices : iterable of int
        1-based indices into `res.exposures` to export.
    scas : iterable of int, optional
        Restrict to a subset of SCAs per exposure. Default: every SCA.
    output : str, optional
        Output CSV filename. If None, auto-generated from the visit_id +
        exposure range.
    show_progress : bool
        Show tqdm progress while streaming.

    Returns
    -------
    str
        Path to the written CSV.
    """
    exposures = [res.select(i) for i in indices]
    if not exposures:
        raise ValueError("No exposures selected — nothing to export.")

    if output is None:
        output = _default_output(res, exposures)

    rows: list = []
    total_files = sum(
        len([s for s in exp.scas if scas is None or s in set(scas)])
        for exp in exposures
    )
    _log(f"Exporting metadata: {len(exposures)} exposure(s), "
         f"~{total_files} SCA file(s)")

    for exp in exposures:
        af_dict = res.stream(exp, scas=scas, show_progress=show_progress)
        try:
            rows.extend(extract_rows(af_dict, exp))
        finally:
            close_streams(af_dict)

    if not rows:
        raise RuntimeError("No metadata extracted — nothing to write.")

    write_csv(rows, output)

    # Consistency check — did we get every exposure the metadata says we should?
    meta_nexp = None
    for r in rows:
        v = r.get('meta.visit.nexposures')
        if v not in (None, ''):
            try:
                meta_nexp = int(v)
                break
            except (TypeError, ValueError):
                pass
    unique_exp_keys = {(r.get('visit_id'), r.get('exposure')) for r in rows}
    if meta_nexp is not None and len(unique_exp_keys) < meta_nexp:
        _log(f"WARNING: meta.visit.nexposures={meta_nexp} but only "
             f"{len(unique_exp_keys)} exposure(s) exported")

    return output


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _resolve_exposure_keys(res: DataResults, spec):
    """Turn --exposures spec (int / 'a,b' / 'a-b' / 'all') into 1-based indices."""
    if spec is None or str(spec).strip().lower() in ('', 'all', '*'):
        return list(range(1, res.n_exposures + 1))
    indices = parse_int_spec(spec)
    for i in indices:
        if i < 1 or i > res.n_exposures:
            raise IndexError(
                f"Exposure index {i} out of range 1..{res.n_exposures} "
                f"(query returned {res.n_exposures} exposure(s))"
            )
    return indices


def _cli():
    import argparse
    from roman_mast import (
        add_list_data_args, list_data_from_args, print_summary,
    )

    p = argparse.ArgumentParser(
        description="Export Roman WFI metadata to a CSV spreadsheet.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Same search flags as roman_fits.py — find the exposures with the standard
--program / --pass / --visit-id / etc., then use --exposures / --scas to
pick which of those to include in the spreadsheet.

Examples:
  # See what's available first (no export yet)
  python export_metadata_csv.py --program 114 --pass 57 --sca-only --list

  # Export every SCA of exposure 1 to a CSV
  python export_metadata_csv.py --program 114 --pass 57 --sca-only \\
      --exposures 1 --output exp1_meta.csv

  # Every exposure in a visit, all SCAs
  python export_metadata_csv.py --visit-id 0011401057001001001 \\
      --sca-only --exposures all

  # Only some SCAs of a range of exposures
  python export_metadata_csv.py --program 114 --pass 57 --sca-only \\
      --exposures 1-3 --scas 1-6

  # Level-1 (uncal) metadata for one visit
  python export_metadata_csv.py --visit-id 0011401057001001001 \\
      --data-level 1 --exposures all
""",
    )

    add_list_data_args(p)

    p.add_argument('--exposures', default='all',
                   help="Which exposure(s) to export (1-based index into the "
                        "listed exposures). '1', '1,3,5', '1-4', or 'all'. "
                        "Default: 'all'.")
    p.add_argument('--scas', default=None,
                   help="Restrict to a subset of SCAs, e.g. '4' / '1-6' / "
                        "'1,3,5'. Default: every SCA the exposure has.")
    p.add_argument('--output', '-o', default=None,
                   help="Output CSV path. Default: auto-generated from "
                        "visit_id + exposure range.")
    p.add_argument('--list', action='store_true',
                   help='List matching exposures and exit without exporting.')
    p.add_argument('--max-rows', type=int, default=50,
                   help='Max exposures to show in the summary (default 50)')

    args = p.parse_args()

    res = list_data_from_args(args)

    if res.n_exposures == 0:
        print_summary(res, max_rows=args.max_rows)
        print("\nNo exposures match — nothing to export.")
        return

    if args.list:
        print_summary(res, max_rows=args.max_rows)
        return

    indices = _resolve_exposure_keys(res, args.exposures)
    scas = parse_int_spec(args.scas) if args.scas is not None else None

    export_csv(res, indices, scas=scas, output=args.output)


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    _cli()
