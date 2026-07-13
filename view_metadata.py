"""Interactive script to stream a Roman WFI exposure and browse its metadata.

This script demonstrates how to:
1. Stream exposure data from MAST
2. Access and explore the roman_datamodels structure
3. Browse metadata hierarchically
4. Export metadata to files
"""

import keyring, keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

import warnings
warnings.filterwarnings('ignore')

import os
import sys
import argparse
import json
import pprint
import roman_datamodels as rdm
from streaming_utils import (
    read_mast_product,
)
from query_utils import add_query_args, prompt_query_params, resolve_query


def get_dict_structure(obj, max_depth=10, current_depth=0, visited=None):
    """Recursively extract the structure of an object as a dictionary.

    Returns a nested dict showing the structure of attributes and their types.
    """
    if visited is None:
        visited = set()

    if current_depth >= max_depth:
        return "<max depth>"

    # Avoid infinite recursion on circular references
    obj_id = id(obj)
    if obj_id in visited:
        return "<circular reference>"
    visited.add(obj_id)

    try:
        # Handle common types
        if isinstance(obj, (str, int, float, bool, type(None))):
            return type(obj).__name__

        if isinstance(obj, (list, tuple)):
            if len(obj) == 0:
                return f"{type(obj).__name__}(empty)"
            return f"{type(obj).__name__}[{len(obj)}]"

        if isinstance(obj, dict):
            return {k: get_dict_structure(v, max_depth, current_depth + 1, visited)
                    for k, v in obj.items()}

        # For objects with __dict__, show their attributes
        if hasattr(obj, '__dict__'):
            result = {}
            for k, v in obj.__dict__.items():
                if not k.startswith('_'):
                    result[k] = get_dict_structure(v, max_depth, current_depth + 1, visited)
            return result if result else type(obj).__name__

        return type(obj).__name__
    except Exception as e:
        return f"<error: {str(e)}>"


def print_metadata_tree(obj, prefix="", max_depth=5, current_depth=0, visited=None):
    """Print metadata in a tree structure."""
    if visited is None:
        visited = set()

    if current_depth >= max_depth:
        print(f"{prefix}...")
        return

    obj_id = id(obj)
    if obj_id in visited:
        print(f"{prefix}<circular>")
        return
    visited.add(obj_id)

    try:
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(v, (dict, list)) and len(str(v)) > 100:
                    print(f"{prefix}{k}: ...")
                else:
                    print(f"{prefix}{k}: {v}")

        elif hasattr(obj, '__dict__'):
            for k, v in sorted(obj.__dict__.items()):
                if not k.startswith('_'):
                    if hasattr(v, '__dict__') or isinstance(v, dict):
                        print(f"{prefix}{k}:")
                        print_metadata_tree(v, prefix + "  ", max_depth, current_depth + 1, visited)
                    elif isinstance(v, (list, tuple)) and len(v) > 0 and hasattr(v[0], '__dict__'):
                        print(f"{prefix}{k}: [{type(v[0]).__name__}] x {len(v)}")
                    else:
                        print(f"{prefix}{k}: {v}")
    except Exception as e:
        print(f"{prefix}<error: {e}>")


def get_value_by_path(obj, path):
    """Get a value from a nested object using dot notation.

    Example: "meta.exposure.start_time"
    """
    parts = path.split('.')
    current = obj

    for part in parts:
        if isinstance(current, dict):
            current = current.get(part)
        else:
            current = getattr(current, part, None)

        if current is None:
            return None

    return current


def browse_metadata_interactive(dm):
    """Interactive metadata browser."""
    print("\n" + "="*80)
    print("METADATA BROWSER")
    print("="*80)
    print("\nCommands:")
    print("  tree                 - Show full metadata tree")
    print("  meta                 - Show top-level metadata groups")
    print("  path <dotted.path>   - Get value at specific path (e.g., 'meta.exposure.start_time')")
    print("  find <keyword>       - Search for keyword in metadata paths")
    print("  export <filename>    - Export metadata to JSON file")
    print("  data                 - Show data array info")
    print("  info                 - Show WCS and basic info")
    print("  help                 - Show this help")
    print("  exit                 - Exit browser")
    print("-"*80 + "\n")

    while True:
        try:
            cmd = input(">>> ").strip()

            if not cmd:
                continue

            if cmd == 'exit':
                break

            elif cmd == 'help':
                print("\nCommands:")
                print("  tree                 - Show full metadata tree")
                print("  meta                 - Show top-level metadata groups")
                print("  path <dotted.path>   - Get value at specific path")
                print("  find <keyword>       - Search for keyword in metadata paths")
                print("  export <filename>    - Export metadata to JSON file")
                print("  data                 - Show data array info")
                print("  info                 - Show WCS and basic info")
                print()

            elif cmd == 'tree':
                print("\n--- Metadata Tree ---")
                print_metadata_tree(dm.meta, max_depth=4)

            elif cmd == 'meta':
                print("\n--- Top-level Metadata Groups ---")
                if hasattr(dm, 'meta'):
                    for attr in sorted(dir(dm.meta)):
                        if not attr.startswith('_'):
                            print(f"  dm.meta.{attr}")

            elif cmd == 'data':
                print(f"\n--- Data Array ---")
                print(f"  Shape: {dm.data.shape}")
                print(f"  Dtype: {dm.data.dtype}")
                print(f"  Size: {dm.data.nbytes / 1e6:.2f} MB")
                print(f"  Min: {dm.data.min():.2f}")
                print(f"  Max: {dm.data.max():.2f}")
                print(f"  Mean: {dm.data.mean():.2f}")

            elif cmd == 'info':
                print(f"\n--- Exposure Info ---")
                paths = [
                    'meta.exposure.start_time',
                    'meta.exposure.end_time',
                    'meta.exposure.effective_exposure_time',
                    'meta.instrument.detector',
                    'meta.instrument.optical_element',
                    'meta.wcsinfo.wcs_type',
                ]
                for path in paths:
                    val = get_value_by_path(dm, path)
                    print(f"  {path}: {val}")

                if hasattr(dm, 'meta') and hasattr(dm.meta, 'wcs'):
                    print(f"\n  WCS Info:")
                    wcs = dm.meta.wcs
                    print(f"    Type: {type(wcs).__name__}")
                    if hasattr(wcs, 'bounding_box'):
                        print(f"    Bounding box: {wcs.bounding_box}")

            elif cmd.startswith('path '):
                path = cmd[5:].strip()
                val = get_value_by_path(dm, path)
                print(f"\n  {path}: {val}")

            elif cmd.startswith('find '):
                keyword = cmd[5:].strip().lower()
                print(f"\n--- Metadata paths containing '{keyword}' ---")

                def find_in_obj(obj, prefix="", matches=None):
                    if matches is None:
                        matches = []

                    if isinstance(obj, dict):
                        for k, v in obj.items():
                            path = f"{prefix}.{k}" if prefix else k
                            if keyword in path.lower() or keyword in str(v).lower():
                                matches.append((path, v))
                            if isinstance(v, (dict, list)):
                                find_in_obj(v, path, matches)

                    elif hasattr(obj, '__dict__'):
                        for k, v in obj.__dict__.items():
                            if not k.startswith('_'):
                                path = f"{prefix}.{k}" if prefix else k
                                if keyword in path.lower():
                                    matches.append((path, v))
                                if hasattr(v, '__dict__') or isinstance(v, dict):
                                    find_in_obj(v, path, matches)

                    return matches

                matches = find_in_obj(dm.meta)
                if matches:
                    for path, val in matches[:20]:  # Limit to 20 results
                        val_str = str(val)[:60]
                        print(f"  meta.{path}: {val_str}")
                    if len(matches) > 20:
                        print(f"  ... and {len(matches) - 20} more")
                else:
                    print("  No matches found")

            elif cmd.startswith('export '):
                filename = cmd[7:].strip()

                # Convert metadata to a serializable dict
                def to_serializable(obj, depth=0, max_depth=5):
                    if depth >= max_depth:
                        return "<max_depth>"

                    if obj is None or isinstance(obj, (str, int, float, bool)):
                        return obj

                    if isinstance(obj, (list, tuple)):
                        return [to_serializable(item, depth+1, max_depth) for item in obj[:10]]

                    if isinstance(obj, dict):
                        return {k: to_serializable(v, depth+1, max_depth) for k, v in obj.items()}

                    if hasattr(obj, '__dict__'):
                        result = {}
                        for k, v in obj.__dict__.items():
                            if not k.startswith('_'):
                                result[k] = to_serializable(v, depth+1, max_depth)
                        return result

                    return str(obj)

                metadata_dict = {
                    'data_shape': dm.data.shape,
                    'data_dtype': str(dm.data.dtype),
                    'metadata': to_serializable(dm.meta)
                }

                with open(filename, 'w') as f:
                    json.dump(metadata_dict, f, indent=2, default=str)

                print(f"  Exported to {filename}")

            else:
                print(f"Unknown command: {cmd}")

        except KeyboardInterrupt:
            print("\nExit")
            break
        except Exception as e:
            print(f"Error: {e}")


def main():
    """Interactive or command-line entry point."""
    parser = argparse.ArgumentParser(
        description='Stream a Roman WFI exposure and browse its metadata interactively.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python view_metadata.py 0012401001001002001
  python view_metadata.py 0012401001001002001 --exp-num 2 --data-level 1 --sca 5
  python view_metadata.py  (interactive mode)
        """,
    )
    add_query_args(parser, visit_wildcard=False, exp_mode='single', sca_mode='single')
    args = parser.parse_args()

    if args.visit_id is None:
        print("Roman WFI Metadata Viewer")
        print("=" * 70)
        params = prompt_query_params(
            visit_wildcard=False,
            exp_mode='single',
            sca_mode='single',
            defaults={'visit_id': '0012401001001002001', 'exp_num': 1, 'sca': 1},
        )
        visit_id = params['visit_id']
        exp_num = int(params['exp_spec']) if params['exp_spec'] else 1
        data_level = params['data_level']
        sca_num = int(params['sca_spec']) if params['sca_spec'] else 1
    else:
        visit_id = args.visit_id
        exp_num = args.exp_num if args.exp_num is not None else 1
        data_level = args.data_level if args.data_level is not None else 2
        sca_num = args.sca if args.sca is not None else 1

    print("\n" + "=" * 70)
    print(f"Streaming {visit_id} exp={exp_num} level={data_level} SCA={sca_num:02d}")
    print("=" * 70)

    # Get filenames via resolve_query so exposure numbers come from filenames
    try:
        q = resolve_query(visit_id, exp_spec=str(exp_num), data_level=data_level,
                          sca_spec=str(sca_num))

        if sca_num not in q.urls:
            print(f"ERROR: SCA {sca_num:02d} not found in results")
            return

        filename = q.urls[sca_num].get(exp_num)
        if filename is None:
            print(f"ERROR: No data available for SCA {sca_num:02d} at exposure {exp_num}")
            return

        # Stream the data
        print(f"\nStreaming SCA {sca_num:02d}...")
        af = read_mast_product(q.missions, filename)

        # Open the data model
        print("Loading data model...")
        dm = rdm.open(af)

        print(f"\nSuccessfully loaded data model")
        print(f"  Shape: {dm.data.shape}")
        print(f"  Dtype: {dm.data.dtype}")

        # Start interactive browser
        browse_metadata_interactive(dm)

        # Cleanup
        af.close()

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
