"""Download and explore a segm_sca ASDF file with fresh S3 URL."""

import sys
import os
import tempfile
import requests
import time

import roman_mast as rm

missions = rm.connect()
filename = "r0011401057001001001_0002_wfi01_f062_segm.asdf"

print(f"Querying for product URL...", file=sys.stderr)
# Query the product list to get fresh S3 URL
res = rm.list_data(
    program=114, pass_=57, observation=1, exposure=2, sca_only=True,
    kinds=['segm_sca'], missions=missions
)

if res.n_products == 0:
    print("No products found", file=sys.stderr)
    sys.exit(1)

wanted = [f for f in res.filenames if f'_wfi01_' in f]
if not wanted:
    print(f"No SCA-01 file found in {res.filenames}", file=sys.stderr)
    sys.exit(1)

filename = wanted[0]
print(f"Found: {filename}", file=sys.stderr)

# Get the S3 URL - this will be fresh
print(f"Getting fresh S3 URL from MAST...", file=sys.stderr)

# Download using the Mastclass's download_file method
with tempfile.TemporaryDirectory() as tmpdir:
    print(f"Downloading to temp dir: {tmpdir}", file=sys.stderr)

    # Download with a fresh signed URL
    try:
        download_path = os.path.join(tmpdir, filename)
        status = missions.download_file(filename, local_path=download_path, verbose=False)
        local_files = [download_path]
        print(f"Retrieved: {local_files}", file=sys.stderr)

        if isinstance(local_files, list) and len(local_files) > 0:
            local_path = local_files[0]
        else:
            local_path = str(local_files)

        print(f"Local path: {local_path}", file=sys.stderr)
        print(f"File exists: {os.path.exists(local_path)}", file=sys.stderr)

        if os.path.exists(local_path):
            import asdf
            print(f"\nOpening ASDF file...\n", file=sys.stderr)

            af = asdf.open(local_path)
            print(f"Type: {type(af)}")

            print("\n" + "=" * 70)
            print("Top-level keys:")
            print("=" * 70)
            for k in sorted(af.tree.keys()):
                print(f"  {k}")

            # Sub-keys
            for meta_key in ('roman', 'meta', 'asdf_library'):
                if meta_key in af.tree:
                    node = af.tree[meta_key]
                    if hasattr(node, 'keys'):
                        print(f"\n[{meta_key}] sub-keys:")
                        for k in sorted(node.keys()):
                            val = node[k]
                            print(f"    {k:<30}  {type(val).__name__}")

            # Try roman_datamodels
            try:
                import roman_datamodels as rdm
                dm = rdm.open(af)
                print(f"\n{'=' * 70}")
                print(f"roman_datamodels type: {type(dm).__name__}")
                print(f"{'=' * 70}\n")

                print("dm attributes:")
                for attr in sorted(dir(dm)):
                    if attr.startswith('_'):
                        continue
                    try:
                        val = getattr(dm, attr)
                        if hasattr(val, 'shape'):
                            print(f"    {attr:<30}  shape={val.shape}  dtype={val.dtype}")
                        elif hasattr(val, '__class__') and 'Node' in type(val).__name__:
                            print(f"    {attr:<30}  (node: {type(val).__name__})")
                    except Exception:
                        pass
                dm.close()
            except Exception as e:
                print(f"\nroman_datamodels failed: {e}")

            af.close()
        else:
            print(f"File does not exist at {local_path}", file=sys.stderr)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

print("\nDone.", file=sys.stderr)
