"""Explore a locally-saved segm_sca ASDF file."""

import sys
import os

# Path to local ASDF file
LOCAL_FILE = "/efs/roman_it_shared/mrizzo/data/r0011401057001001001_0002_wfi01_f062_segm.asdf"

if not os.path.exists(LOCAL_FILE):
    print(f"File not found: {LOCAL_FILE}", file=sys.stderr)
    sys.exit(1)

print(f"Opening {LOCAL_FILE}...", file=sys.stderr)

import asdf

af = asdf.open(LOCAL_FILE)
print(f"Type returned: {type(af)}")

print()
print("=" * 70)
print("ASDF tree keys (top level):")
print("=" * 70)
for k in sorted(af.tree.keys()):
    print(f"  {k}")

# If there's a roman_meta / meta, show one more level
for meta_key in ('roman', 'meta', 'asdf_library'):
    if meta_key in af.tree:
        node = af.tree[meta_key]
        if hasattr(node, 'keys'):
            print(f"\n[{meta_key}] sub-keys:")
            for k in sorted(node.keys()):
                val = node[k]
                print(f"    {k:<30}  {type(val).__name__}")

# Try opening as a roman datamodel for structured access
try:
    import roman_datamodels as rdm
    dm = rdm.open(af)
    print(f"\n{'=' * 70}")
    print(f"roman_datamodels type: {type(dm).__name__}")
    print(f"{'=' * 70}")

    # Show top-level attributes that have data
    print("\ndm attributes (with shapes/types):")
    for attr in sorted(dir(dm)):
        if attr.startswith('_'):
            continue
        try:
            val = getattr(dm, attr)
            if hasattr(val, 'shape'):
                print(f"    {attr:<30}  shape={val.shape}  dtype={val.dtype}")
            elif hasattr(val, '__class__') and 'Node' in type(val).__name__:
                print(f"    {attr:<30}  (node: {type(val).__name__})")
        except Exception as e:
            pass
    dm.close()
except Exception as e:
    print(f"\nroman_datamodels.open failed: {e}")

af.close()
print("\nDone.")
