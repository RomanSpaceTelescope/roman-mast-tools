"""Check whether segm unique IDs match cat_sca row count for SCA 01."""

import sys, os, tempfile, warnings
warnings.filterwarnings('ignore')

import roman_mast as rm
import numpy as np

missions = rm.connect()

TARGET = dict(program=114, pass_=57, observation=1, exposure=2, sca_only=True)
SCA_PAT = '_wfi01_'

print("Querying segm_sca ...", flush=True)
res_segm = rm.list_data(**TARGET, kinds=['segm_sca'], missions=missions)
segm_file = next((f for f in res_segm.filenames if SCA_PAT in f), None)
print(f"  segm: {segm_file}")

print("Querying cat_sca ...", flush=True)
res_cat = rm.list_data(**TARGET, kinds=['cat_sca'], missions=missions)
cat_file = next((f for f in res_cat.filenames if SCA_PAT in f), None)
print(f"  cat : {cat_file}")

if not segm_file or not cat_file:
    print("ERROR: could not find both files")
    sys.exit(1)

with tempfile.TemporaryDirectory() as tmpdir:
    # Download segm
    print(f"\nDownloading segm ...", flush=True)
    segm_path = os.path.join(tmpdir, segm_file)
    missions.download_file(segm_file, local_path=segm_path, verbose=False)

    import asdf
    af = asdf.open(segm_path)
    data = af.tree['roman']['data']
    unique_ids = np.unique(data)
    n_nonzero = int((unique_ids > 0).sum())
    print(f"  unique IDs in segm (incl. 0=background): {len(unique_ids)}")
    print(f"  non-background segment IDs: {n_nonzero}")
    af.close()

    # Download cat
    print(f"\nDownloading cat ...", flush=True)
    cat_path = os.path.join(tmpdir, cat_file)
    missions.download_file(cat_file, local_path=cat_path, verbose=False)

    import pandas as pd
    df = pd.read_parquet(cat_path)
    print(f"  cat rows: {len(df)}")

print(f"\nMatch (non-zero segm IDs == cat rows): {n_nonzero == len(df)}")
print(f"  segm non-zero IDs: {n_nonzero}")
print(f"  cat rows:          {len(df)}")
print(f"  difference:        {abs(n_nonzero - len(df))}")
