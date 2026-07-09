"""Export all metadata mnemonics from Roman WFI exposures as a CSV table.

Creates a CSV where:
- Each row is an SCA (or exposure/SCA combination)
- Each column is a metadata mnemonic (field key)
- Values are the metadata values
"""

import keyring, keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

import warnings
warnings.filterwarnings('ignore')

import os
import argparse
import csv
import io
import json
from collections import defaultdict
import roman_datamodels as rdm
from streaming_utils import (
    get_MAST_token,
    group_file_urls,
    stream_to_buffer,
    close_buffer_streams,
)
from query_utils import add_query_args, prompt_query_params, resolve_query


def flatten_metadata(obj, prefix="", visited=None, max_depth=10, current_depth=0):
    """Recursively flatten nested metadata into dot-notation keys.

    Returns a dict of {key: value} where keys use dot notation.
    Skips private attributes and deeply nested objects.
    Handles DNode objects from roman_datamodels.
    """
    if visited is None:
        visited = set()

    if current_depth >= max_depth:
        return {}

    result = {}
    obj_id = id(obj)

    if obj_id in visited:
        return {}
    visited.add(obj_id)

    try:
        # Handle dict-like objects (including DNode which supports __getitem__)
        if isinstance(obj, dict) or (hasattr(obj, '__getitem__') and hasattr(obj, 'keys')):
            try:
                items = obj.items()
            except:
                try:
                    items = [(k, obj[k]) for k in obj.keys()]
                except:
                    return {}

            for k, v in items:
                key = f"{prefix}.{k}" if prefix else k

                # Skip None values
                if v is None:
                    continue

                # Skip array data but include metadata
                if isinstance(v, (list, tuple)):
                    # For lists of primitives, convert to string
                    if len(v) == 0:
                        pass
                    elif all(isinstance(x, (str, int, float, bool, type(None))) for x in v):
                        result[key] = str(v)[:200]  # Limit length
                    elif len(v) < 20:
                        # Recurse into smaller lists
                        nested = flatten_metadata(v, key, visited, max_depth, current_depth + 1)
                        result.update(nested)
                    # Skip large lists
                    continue

                # Handle scalar types
                if isinstance(v, (str, int, float, bool)):
                    result[key] = v

                # Handle Time objects and other special types
                elif hasattr(v, 'iso') or hasattr(v, 'datetime'):
                    # Likely an astropy Time object
                    result[key] = str(v)

                # Handle dict-like and nested objects
                elif hasattr(v, '__getitem__') or hasattr(v, '__dict__'):
                    nested = flatten_metadata(v, key, visited, max_depth, current_depth + 1)
                    result.update(nested)

        elif hasattr(obj, '__dict__'):
            for k, v in obj.__dict__.items():
                if k.startswith('_'):
                    continue

                key = f"{prefix}.{k}" if prefix else k

                if v is None:
                    continue

                if isinstance(v, (str, int, float, bool)):
                    result[key] = v
                elif isinstance(v, (list, tuple)):
                    if len(v) == 0:
                        pass
                    elif all(isinstance(x, (str, int, float, bool, type(None))) for x in v):
                        result[key] = str(v)[:200]
                    elif len(v) < 20:
                        nested = flatten_metadata(v, key, visited, max_depth, current_depth + 1)
                        result.update(nested)
                    continue

                elif hasattr(v, 'iso') or hasattr(v, 'datetime'):
                    result[key] = str(v)

                elif hasattr(v, '__getitem__') or hasattr(v, '__dict__'):
                    nested = flatten_metadata(v, key, visited, max_depth, current_depth + 1)
                    result.update(nested)

    except Exception as e:
        pass

    return result


def extract_metadata_from_sca(buf, sca_num):
    """Open an SCA buffer and extract all metadata as a flat dict."""
    try:
        import numpy as np

        buf.seek(0)
        dm = rdm.open(buf)

        # Load data and compute stats (handle NaN values)
        data_array = np.asarray(dm.data[...])

        # Extract visit ID from filename if available
        visit_id = None
        if hasattr(dm, 'meta') and hasattr(dm.meta, 'filename'):
            filename = dm.meta.filename
            if filename:
                # Roman filenames start with 'r' followed by visit ID
                # Format: r<visit_id>_<exp>_<sca>_<filter>_<level>
                try:
                    visit_id = filename.split('_')[0].replace('r', '')
                except:
                    pass

        # Get basic info
        row = {
            'visit': visit_id or '',
            'sca': sca_num,
            'data_shape': str(dm.data.shape),
            'data_dtype': str(dm.data.dtype),
            'data_min': float(np.nanmin(data_array)),
            'data_max': float(np.nanmax(data_array)),
            'data_mean': float(np.nanmean(data_array)),
            'data_valid_pixels': int(np.isfinite(data_array).sum()),
            'data_nan_pixels': int(np.isnan(data_array).sum()),
        }

        # Flatten all metadata
        if hasattr(dm, 'meta'):
            meta_flat = flatten_metadata(dm.meta, prefix='meta')
            row.update(meta_flat)

        return row

    except Exception as e:
        print(f"    ERROR reading SCA {sca_num}: {e}")
        return {'sca': sca_num, 'error': str(e)}


def export_metadata_csv(visit_id, exp_num=None, exp_num_start=None, exp_num_end=None, data_level=2, output_file=None, scas=None):
    """Export metadata from all SCAs to a CSV file.

    Args:
        visit_id (str): Visit ID to query. Supports wildcards (e.g., '001240100100100*').
        exp_num (int, optional): Single exposure number to query. Mutually exclusive with exp_num_start/end.
        exp_num_start (int, optional): Starting exposure number (inclusive)
        exp_num_end (int, optional): Ending exposure number (inclusive)
                                    If all exposure params are None, all available exposures are exported.
        data_level (int): Data level (1=uncal, 2=cal)
        output_file (str): Output CSV filename. If None, auto-generates.
        scas (list): List of SCA numbers to export. If None, exports all available.

    Returns:
        str: Path to output CSV file
    """

    # Build exp_spec string for resolve_query
    if exp_num is not None:
        exp_spec = str(exp_num)
    elif exp_num_start is not None and exp_num_end is not None:
        exp_spec = f"{exp_num_start}-{exp_num_end}"
    else:
        exp_spec = None  # all available

    sca_spec = ','.join(str(s) for s in scas) if scas else None

    q = resolve_query(visit_id, exp_spec=exp_spec, data_level=data_level, sca_spec=sca_spec)

    mast_token = get_MAST_token()
    urls = q.urls
    exp_nums = q.exp_nums
    scas = q.scas

    if not scas:
        raise ValueError("No valid SCAs to export")

    # Stream and extract metadata from each SCA and exposure
    all_rows = []
    buffer_dict = {}
    meta_visit_nexposures = None  # Will be set from first streamed file

    total_to_stream = len(scas) * len(exp_nums)
    print(f"\nStreaming {total_to_stream} file(s) ({len(scas)} SCAs × {len(exp_nums)} exposures)...")

    stream_count = 0
    for exp_num in exp_nums:
        print(f"\n  Exposure {exp_num}:")

        for sca_num in scas:
            try:
                url = urls[sca_num].get(exp_num)
                if url is None:
                    print(f"    SCA {sca_num:02d}: No data available")
                    continue

                stream_count += 1
                print(f"    SCA {sca_num:02d}: Streaming...", end=" ", flush=True)
                buf = stream_to_buffer(url, mast_token=mast_token, show_progress=False)

                print("Extracting metadata...", end=" ", flush=True)
                row = extract_metadata_from_sca(buf, sca_num)
                # Add exposure number to each row
                row['exposure'] = exp_num
                all_rows.append(row)
                buffer_dict[f"exp{exp_num}_sca{sca_num}"] = buf
                print("✓")

                # Check meta.visit.nexposures on first extraction
                if meta_visit_nexposures is None and 'meta.visit.nexposures' in row:
                    try:
                        meta_visit_nexposures = int(row['meta.visit.nexposures'])
                    except (ValueError, TypeError):
                        pass

            except Exception as e:
                print(f"ERROR: {e}")

    # Check if actual exposures is less than meta.visit.nexposures
    if meta_visit_nexposures is not None and len(exp_nums) < meta_visit_nexposures:
        print("\n" + "!"*80)
        print("⚠️  WARNING: INCOMPLETE EXPOSURE SET")
        print("!"*80)
        print(f"  meta.visit.nexposures: {meta_visit_nexposures}")
        print(f"  Actual exposures available: {len(exp_nums)}")
        print(f"  Missing: {meta_visit_nexposures - len(exp_nums)} exposure(s)")
        print("!"*80 + "\n")

    # Combine all rows and collect all unique keys
    all_keys = set()
    for row in all_rows:
        all_keys.update(row.keys())

    all_keys = sorted(all_keys)

    # Reorder keys with user-specified priority order
    priority_keys = [
        'visit', 'exposure', 'sca',
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
        'data_shape', 'data_dtype', 'data_min', 'data_max', 'data_mean', 'data_nan_pixels', 'data_valid_pixels',
        'meta.source_catalog.tweakreg_catalog_name',

    ]

    ordered_keys = [k for k in priority_keys if k in all_keys]
    remaining_keys = sorted([k for k in all_keys if k not in ordered_keys])
    ordered_keys.extend(remaining_keys)

    # Write CSV
    if output_file is None:
        if len(exp_nums) == 1:
            output_file = f"metadata_{visit_id}_exp{exp_nums[0]:02d}_level{data_level}.csv"
        else:
            output_file = f"metadata_{visit_id}_exp{exp_nums[0]:02d}-{exp_nums[-1]:02d}_level{data_level}.csv"

    print(f"\nWriting CSV: {output_file}")
    print(f"  Rows (SCA/Exposure pairs): {len(all_rows)}")
    print(f"  Columns (mnemonics): {len(ordered_keys)}")

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=ordered_keys, restval='')
        writer.writeheader()

        for row in all_rows:
            # Fill missing values
            for key in ordered_keys:
                if key not in row:
                    row[key] = ''
            writer.writerow(row)

    print(f"\n✓ Exported to: {output_file}")

    # Print summary
    print(f"\nMetadata Summary:")
    if len(all_rows) == 0:
        print(f"  ERROR: No data rows extracted!")
        return None

    # Count unique visits and exposures in results
    unique_visits = set()
    unique_exps = set()
    for row in all_rows:
        # Extract visit from filename if available
        if 'meta.filename' in row:
            filename = row['meta.filename']
            if filename:
                # Roman filenames start with visit ID
                visit_part = filename.split('_')[0].replace('r', '')
                unique_visits.add(visit_part)
        unique_exps.add(row.get('exposure', '?'))

    print(f"  Total rows (exposure/SCA pairs): {len(all_rows)}")
    print(f"  Unique visits: {len(unique_visits)} ({', '.join(sorted(unique_visits)[:3])}{'...' if len(unique_visits) > 3 else ''})")
    print(f"  Total exposures processed: {len(exp_nums)}")
    print(f"  SCAs per exposure: {len(scas)}")
    print(f"  Total mnemonics: {len(ordered_keys)}")
    print(f"\nFirst 20 mnemonics:")
    for i, key in enumerate(ordered_keys[:20]):
        print(f"  {i+1:2d}. {key}")

    if len(ordered_keys) > 20:
        print(f"  ... and {len(ordered_keys) - 20} more")

    # Cleanup buffers
    close_buffer_streams(buffer_dict)

    return output_file


def main():
    """Interactive or command-line entry point."""
    parser = argparse.ArgumentParser(
        description='Export Roman WFI metadata to a CSV file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python export_metadata_csv.py 0012401001001002001
  python export_metadata_csv.py 0012401001001002001 --exp-range 1-3 --data-level 2
  python export_metadata_csv.py '001240100100100*' --scas 1,2,3 --output metadata.csv
  python export_metadata_csv.py  (interactive mode)
        """,
    )
    add_query_args(parser, visit_wildcard=False, exp_mode='flexible', sca_mode='all')
    parser.add_argument('--output', default=None,
                        help='Output CSV filename (default: auto-generated)')
    args = parser.parse_args()

    if args.visit_id is None:
        print("Roman WFI Metadata CSV Exporter")
        print("=" * 70)
        params = prompt_query_params(
            visit_wildcard=False,
            exp_mode='flexible',
            sca_mode='all',
        )
        visit_id = params['visit_id']
        exp_spec = params['exp_spec']
        data_level = params['data_level']
        sca_spec = params['sca_spec']
        output_file = input("Enter output filename (blank = auto-generated): ").strip() or None
    else:
        visit_id = args.visit_id
        exp_spec = getattr(args, 'exp_range', None) or (str(args.exp_num) if getattr(args, 'exp_num', None) else None)
        data_level = args.data_level if args.data_level is not None else 2
        sca_spec = args.scas
        output_file = args.output

    # Translate exp_spec into the legacy kwargs expected by export_metadata_csv()
    from query_utils import parse_exp_spec
    exp_list = parse_exp_spec(exp_spec)
    if exp_list is None:
        exp_num = exp_num_start = exp_num_end = None
    elif len(exp_list) == 1:
        exp_num = exp_list[0]
        exp_num_start = exp_num_end = None
    else:
        exp_num = None
        exp_num_start = exp_list[0]
        exp_num_end = exp_list[-1]

    from query_utils import parse_sca_spec
    scas = parse_sca_spec(sca_spec)

    try:
        csv_file = export_metadata_csv(
            visit_id=visit_id,
            exp_num=exp_num,
            exp_num_start=exp_num_start,
            exp_num_end=exp_num_end,
            data_level=data_level,
            output_file=output_file,
            scas=scas,
        )
        print(f"\nSuccess! File saved to: {csv_file}")

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
