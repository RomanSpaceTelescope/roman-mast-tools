"""Sandbox: explore wcs / segm_sca / cat_sca products for one exposure.

Usage:
    conda run -n roman-mast-tools python explore_products.py

Streams a single SCA of each product kind and dumps its structure so we can
see what fields / arrays are available before wiring anything into the stack.
"""

import sys
import os
import warnings
warnings.filterwarnings('ignore')

import roman_mast as rm

# ── narrow query: one exposure, one SCA ──────────────────────────────────────
TARGET = dict(program=114, pass_=57, observation=1, exposure=2, sca_only=True)
SCA    = 1   # just one SCA to keep streaming fast

kinds_to_explore = ['segm_sca']

print(f"Connecting to MAST ...", file=sys.stderr)
missions = rm.connect()

# ── stream each kind independently ───────────────────────────────────────────
for kind in kinds_to_explore:
    print()
    print("=" * 70)
    print(f"  KIND: {kind}")
    print("=" * 70)

    res = rm.list_data(**TARGET, kinds=[kind], missions=missions)
    if res.n_products == 0:
        print("  (no products found)")
        continue

    # Find the filename for the requested SCA.
    wanted = [f for f in res.filenames if f'_wfi{SCA:02d}_' in f]
    if not wanted:
        print(f"  No filename found for SCA {SCA:02d} in {res.filenames[:3]}")
        continue
    filename = wanted[0]
    print(f"  Filename : {filename}")

    # ── ASDF products ────────────────────────────────────────────────────────
    if filename.endswith('.asdf'):
        print(f"  Streaming via missions.read_product (with fresh S3 URL) ...")

        try:
            af = missions.read_product(filename)
            print(f"  Type returned: {type(af)}")

            print()
            print("  ASDF tree keys (top level):")
            for k in sorted(af.tree.keys()):
                print(f"    {k}")

            # If there's a roman_meta / meta, show one more level
            for meta_key in ('roman', 'meta', 'asdf_library'):
                if meta_key in af.tree:
                    node = af.tree[meta_key]
                    if hasattr(node, 'keys'):
                        print(f"\n  [{meta_key}] sub-keys:")
                        for k in sorted(node.keys()):
                            val = node[k]
                            print(f"    {k:<30}  {type(val).__name__}")

            # Try opening as a roman datamodel for structured access
            try:
                import roman_datamodels as rdm
                dm = rdm.open(af)
                print(f"\n  roman_datamodels type: {type(dm).__name__}")

                # Show top-level attributes that have data
                print("  dm attributes:")
                for attr in dir(dm):
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
                print(f"\n  roman_datamodels.open failed: {e}")

            af.close()
        except Exception as e:
            print(f"  download_products failed: {e}")

    # ── Parquet products ─────────────────────────────────────────────────────
    elif filename.endswith('.parquet'):
        # read_product only supports .fits/.asdf — download via get_product_list
        # + download_products, or hit the S3 URL directly via requests.
        print(f"  read_product unsupported for parquet; trying download_products ...")
        import tempfile, os, io
        import pandas as pd

        # missions.download_products expects a products Table.
        # Re-use the filtered table from list_data if it has a 'filename' column.
        from astropy.table import Table
        prod_table = Table({'filename': [filename]})
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                manifest = missions.download_products(prod_table, download_dir=tmpdir, flat=True)
                print(f"  download manifest:\n  {manifest}")
                # Find the downloaded file
                local_path = None
                for row in manifest:
                    if str(row.get('Local Path', '')).endswith('.parquet'):
                        local_path = row['Local Path']
                        break
                if local_path and os.path.exists(local_path):
                    df = pd.read_parquet(local_path)
                    result = df
                else:
                    print("  Could not locate downloaded parquet file in manifest")
                    result = None
        except Exception as e:
            print(f"  download_products failed: {e}")
            result = None

        if isinstance(result, pd.DataFrame):
            df = result
            print(f"\n  rows   : {len(df)}")
            print(f"  columns: {len(df.columns)}")
            print()
            print("  Columns + dtypes:")
            for col in df.columns:
                print(f"    {col:<40}  {df[col].dtype}")
            print()
            print("  First row:")
            print(df.iloc[0].to_string())

    else:
        print(f"  Unknown extension — raw repr:")
        result = missions.read_product(filename)
        print(f"  type: {type(result)}")
        print(f"  repr: {repr(result)[:300]}")

print()
print("Done.")
