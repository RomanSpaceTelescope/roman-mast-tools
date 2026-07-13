"""Detailed exploration of segm_sca ASDF structure and contents."""

import sys
import os
import tempfile

import roman_mast as rm
import asdf
import roman_datamodels as rdm
import numpy as np

missions = rm.connect()
filename = "r0011401057001001001_0002_wfi01_f062_segm.asdf"

print(f"Querying for {filename}...", file=sys.stderr)
res = rm.list_data(
    program=114, pass_=57, observation=1, exposure=2, sca_only=True,
    kinds=['segm_sca'], missions=missions
)

wanted = [f for f in res.filenames if f'_wfi01_' in f]
if not wanted:
    print("No SCA-01 file found", file=sys.stderr)
    sys.exit(1)

filename = wanted[0]
print(f"Found: {filename}\n", file=sys.stderr)

with tempfile.TemporaryDirectory() as tmpdir:
    # Download
    download_path = os.path.join(tmpdir, filename)
    status = missions.download_file(filename, local_path=download_path, verbose=False)

    if not os.path.exists(download_path):
        print(f"Download failed", file=sys.stderr)
        sys.exit(1)

    # Open both ASDF and datamodel
    af = asdf.open(download_path)
    dm = rdm.open(af)

    print("=" * 80)
    print(f"File: {filename}")
    print(f"Data model type: {type(dm).__name__}")
    print("=" * 80)

    # Data arrays
    print("\nData Arrays:")
    print("-" * 80)
    print(f"  data              : shape={dm.data.shape}, dtype={dm.data.dtype}")
    print(f"  detection_image   : shape={dm.detection_image.shape}, dtype={dm.detection_image.dtype}")

    # Stats on the arrays
    print(f"\n  data stats:")
    print(f"    min={dm.data.min()}, max={dm.data.max()}")
    print(f"    unique values: {len(np.unique(dm.data))}")

    print(f"\n  detection_image stats:")
    print(f"    min={dm.detection_image.min():.3f}, max={dm.detection_image.max():.3f}")

    # Metadata
    print("\nMetadata Structure:")
    print("-" * 80)

    def print_dict(d, prefix="  ", depth=0, max_depth=3):
        """Recursively print dict/node structure."""
        if depth >= max_depth:
            return
        for key in sorted(d.keys() if hasattr(d, 'keys') else []):
            try:
                val = d[key]
                if hasattr(val, 'keys'):  # dict-like
                    print(f"{prefix}{key}:")
                    print_dict(val, prefix + "  ", depth + 1, max_depth)
                elif isinstance(val, (list, tuple)):
                    print(f"{prefix}{key}: <{type(val).__name__}, len={len(val)}>")
                elif isinstance(val, np.ndarray):
                    print(f"{prefix}{key}: <ndarray shape={val.shape}, dtype={val.dtype}>")
                else:
                    val_str = str(val)[:60]
                    print(f"{prefix}{key}: {val_str}")
            except Exception as e:
                print(f"{prefix}{key}: <error: {e}>")

    print_dict(dm.meta, depth=0, max_depth=4)

    print("\n" + "=" * 80)
    print("Done.")
    print("=" * 80)

    dm.close()
    af.close()
