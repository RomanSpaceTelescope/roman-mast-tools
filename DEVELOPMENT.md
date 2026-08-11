# Development Workflow for roman_view_sca.py

## Quick Start: Work Locally with Cached Data

To avoid re-streaming from S3 every time you change code, download a sample ASDF file once and work locally.

### 1. Download and Cache a Sample File

```bash
conda run -n roman-mast-tools python cache_sca.py \
    s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
    r0003201001001001004_0001_wfi11_f106_cal.asdf
```

This downloads the file to `./cache/` (default). To use a different cache directory:

```bash
conda run -n roman-mast-tools python cache_sca.py \
    s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
    r0003201001001001004_0001_wfi11_f106_cal.asdf \
    --cache-dir ~/roman_dev_cache
```

The script will print the command you should use next.

### 2. Develop Locally

Once cached, iterate quickly on code changes using the local copy:

```bash
# Display with matplotlib (fast iteration)
conda run -n roman-mast-tools python roman_view_sca.py \
    ./cache/ r0003201001001001004_0001_wfi11_f106_cal.asdf \
    --display mpl

# With photometry overlay
conda run -n roman-mast-tools python roman_view_sca.py \
    ./cache/ r0003201001001001004_0001_wfi11_f106_cal.asdf \
    --display mpl --phot --phot-out phot.csv

# With channel dividers
conda run -n roman-mast-tools python roman_view_sca.py \
    ./cache/ r0003201001001001004_0001_wfi11_f106_cal.asdf \
    --display mpl --channels --grid
```

### 3. Cache Multiple SCAs (Optional)

For broader testing, cache files from different SCAs or filters:

```bash
# SCA 11, F106
conda run -n roman-mast-tools python cache_sca.py \
    s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
    r0003201001001001004_0001_wfi11_f106_cal.asdf --cache-dir ~/roman_dev_cache

# SCA 14, F212
conda run -n roman-mast-tools python cache_sca.py \
    s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
    r0003201001001001004_0001_wfi14_f212_cal.asdf --cache-dir ~/roman_dev_cache
```

Then switch between them during development:

```bash
conda run -n roman-mast-tools python roman_view_sca.py \
    ~/roman_dev_cache/ r0003201001001001004_0001_wfi14_f212_cal.asdf --display mpl
```

## How It Works

- **`stream_sca(uri, filename)`** now detects local vs. S3 paths automatically:
  - If `uri` starts with `/` (local path): opens the file directly with `open()`
  - Otherwise: treats it as an S3 URI and uses `s3fs.S3FileSystem(anon=True)`

- **`cache_sca.py`** is a lightweight wrapper around `s3fs.get()` for one-time downloads

- No changes to photometry, display, or region logic—all improvements apply to both local and S3 use

## Testing End-to-End

Once you're satisfied with local changes, verify against live S3 (reference workflow):

```bash
# Live S3 (as in the existing tutorial)
conda run -n roman-mast-tools python roman_view_sca.py \
    s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
    r0003201001001001004_0001_wfi11_f106_cal.asdf \
    --display mpl --phot
```

This confirms your changes work with real S3 streaming as well.
