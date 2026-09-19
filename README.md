# roman-mast-tools

Tools for querying, streaming, and analyzing Roman Space Telescope WFI data
from MAST — plus telemetry queries against the MAST Engineering Database.

Everything streams directly from S3/MAST into memory; nothing has to be
downloaded first (though a local cache is supported for fast iteration).

---

## Installation

Create and activate the conda environment:

```bash
mamba env create -f environment.yml
mamba activate roman-mast-tools
```

If the environment already exists and you need to update it:

```bash
mamba env update -f environment.yml --prune
```

Verify the install:

```bash
python -c "
import roman_datamodels, romancal, astroquery, fsspec, matplotlib
from astroquery.mast import MastMissions
print('roman_datamodels', roman_datamodels.__version__)
print('romancal        ', romancal.__version__)
print('astroquery      ', astroquery.__version__)
print('OK')
"
```

All commands below assume the `roman-mast-tools` environment is active.

This installs the following command-line tools:

| Command | Purpose |
| --- | --- |
| `roman-mast` | Query and list Roman observations from MAST |
| `roman-fits` | Stream WFI SCAs to FITS files or a DS9 mosaic |
| `roman-metadata` | Export per-SCA ASDF metadata to CSV |
| `roman-view-sca` | Stream and visualize a single SCA (DS9 or matplotlib) |
| `roman-phot` | Batch aperture photometry across all SCAs of an exposure |
| `roman-color` | Build an RGB composite of one SCA from three filters and display it in DS9 |
| `roman-telem` | Query Roman telemetry mnemonics from the MAST Engineering DB |
| `roman-telem-plot` | Plot telemetry CSV/Parquet with optional grouping |

---

## Authentication

MAST authentication is only required for **private programs** and **all
telemetry queries**. Public observation data (including the tutorial
dataset) is accessed anonymously via S3.

### MAST observations (`roman-mast`, `roman-fits`, `roman-view-sca`, `roman-phot`, `roman-metadata`)

```bash
export MAST_API_TOKEN=your_token_here
```

### Telemetry (`roman-telem`)

The MAST Engineering Database requires a token:

```bash
export MAST_API_TOKEN=your_token_here
```

### Getting a token

1. Visit the [MAST Portal](https://auth.mast.stsci.edu/)
2. Sign in with your STScI account
3. Go to **Settings → API Tokens**
4. Generate a new token and save it securely

---

## Quick Start

```bash
# 1. List available exposures
roman-mast --program 1039

# 2. Stream exposure 1 to FITS files
# First you need to run "ds9&" in your terminal
roman-fits --program 1039 --execution-plan 2 --pass 4 --observation 7 \
    --exposures 1 --to ds9

# 3. View a single SCA using matplotlib
roman-view-sca --display mpl \
    --program 1039 --execution-plan 2 --pass 4 --observation 7 --exposure 1 --sca 11
```

---

## `roman-mast` — Query and Stream from MAST

`roman_mast.py` is the foundation of the stack. Every other tool builds on
its streaming primitives.

### Basic usage

```bash
# List all exposures for a program/pass
roman-mast --program 1039 --execution-plan 2 --pass 4 --observation 7

# List only the per-SCA L2 cal files for a specific detector
roman-mast --program 1039 --execution-plan 2 --pass 4 --observation 7 --detector wfi04 --sca-only --list

# See every supported product kind (cal, uncal, cat, cat_sca, coadd, gw, ...)
roman-mast --list-kinds

# Fetch a specific set of product kinds
roman-mast --program 1039 --execution-plan 2 --pass 4 --observation 7 --kinds cal,cat_sca --list
```

### Python API

```python
from roman_mast import list_data, print_summary

# All filters are optional
res = list_data(program=1039, execution_plan=2, pass_=4, observation=7, sca_only=True)
print_summary(res)

# Stream one exposure into memory — key is a 1-based index or (visit_id, exposure) tuple
af_dict = res.stream(1)                       # 1-based index from print_summary
# af_dict = res.stream(('0103900004001001007', 1))   # or tuple form

# One-liner: stream → FITS files on disk
res.to_fits(1, out_dir='output/')

# One-liner: stream → DS9 mosaic (DS9 must be running)
res.to_ds9(1)
```

`res.stream(...)` returns `{sca: AsdfFile}`. `res.to_fits(...)` and
`res.to_ds9(...)` accept the same key. See `roman_fits.py` for details.

### Filter reference

| Filter | Type | Description |
| --- | --- | --- |
| `program` | int | APT program ID (e.g. `114`) |
| `pass_` | int | Pass within the execution plan (kwarg is `pass_` because `pass` is a Python keyword; on the CLI it's `--pass`) |
| `execution_plan` | int | Execution plan within the program |
| `segment` | int | Segment within the pass |
| `observation` | int | Observation within the segment |
| `visit` | int | Visit within the observation |
| `detector` | str/int | `'WFI04'`, `'wfi04'`, or `4` |
| `visit_id` | str | Full 19-digit ID or wildcard, e.g. `'0103900004*'` |
| `exposure` | int | Matches last 4 digits of `observation_id` |
| `optical_element` | str | e.g. `'F062'`, `'F129'` |
| `exposure_type` | str | e.g. `'WFI_IMAGE'`, `'WFI_DARK'` |
| `product_type` | str | `'l2'` (per-SCA) or `'p_visit_coadd'` (mosaic tiles) |
| `sca_only` | bool | Shortcut for `product_type='l2'`; drops mosaic tiles |
| `data_level` | 1 / 2 / `'gw'` / `None` | 1 → `_uncal`, 2 → `_cal`, `'gw'` → `_gw`, `None` → all kinds |
| `kinds` | list or `'all'` | Explicit product-kind list (see `--list-kinds`). Overrides `data_level`. |

### Useful CLI options

- `--list` — print matching exposures and exit
- `--kinds cal,cat_sca` / `--kinds all` — request specific product kinds
- `--list-kinds` — show every supported product kind
- `--enumerate-products` — hit MAST for the authoritative product list (slow; usually unnecessary)
- `--server URL` — override MAST server URL
- `--quiet`, `--color {auto,always,never}`, `--no-color`

---

## `roman-fits` — Stream WFI SCAs to FITS Files or DS9

`roman_fits.py` builds on `roman_mast` to stream all 18 SCAs of an exposure
straight into either FITS files (with SIP-approximated WCS headers) or a
DS9 WCS mosaic — plus optional L4 catalog overlays.

### Command line

Selection is two-step: `--program` / `--pass` / `--visit-id` / etc. find
the matching **exposures**, then `--exposures` picks which of those to
output.

```bash
# See what's available (no output yet)
roman-fits --program 1039 --execution-plan 2 --pass 4 --observation 7 --sca-only --list

# Write one exposure to disk as 18 FITS files under /tmp/wfi/v..._exp01/
roman-fits --program 1039 --execution-plan 2 --pass 4 --observation 7 --sca-only \
    --exposures 1 --to fits --out-dir /tmp/wfi

# RICE compressed
roman-fits --program 1039 --execution-plan 2 --pass 4 --observation 7 --sca-only \
    --exposures 1 --to fits --compress

# Multiple exposures — one v..._expNN/ folder each under --out-dir
roman-fits --program 1039 --execution-plan 2 --pass 4 --observation 7 --sca-only \
    --exposures 1-4 --to fits --out-dir /tmp/wfi

# Stream to DS9 (DS9 + pyds9 required, `ds9 &` running). Catalog parquet
# files + a combined .reg land in <cwd>/v..._exp01/catalog/
roman-fits --program 1039 --execution-plan 2 --pass 4 --observation 7 --sca-only \
    --exposures 1 --to ds9

# Only a subset of SCAs
roman-fits --program 1039 --execution-plan 2 --pass 4 --observation 7 --sca-only \
    --exposures 1 --scas 1-6 --to fits

# Coadd tiles (product_type=p_visit_coadd) instead of per-SCA exposures
roman-fits --program 1039 --execution-plan 2 --pass 4 --observation 7 --coadd --exposures 1 --to fits
```

### Output layout

Each run drops **one folder per exposure** under `--out-dir` (default cwd):

```
<out-dir>/v{visit_id}_exp{NN}/
    metadata_{visit_id}_exp{NN}.csv    # unless --no-metadata
    catalog/                            # ds9 mode, unless --no-catalog
        r..._wfi{NN}_..._cat.parquet    # one per SCA that has a catalog
        catalog_{visit_id}_exp{NN}.reg  # DS9 regions from all catalogs
    sca_{NN}.fits                       # fits mode
```

### Useful options

- `--to {fits,ds9}` — output mode (default `fits`)
- `--exposures 1` / `1,3,5` / `1-4` / `all` — which matching exposures to output
- `--scas 4` / `1-6` / `1,3,5` — restrict to a subset of SCAs
- `--compress` — write RICE-tiled `.fits.fz`
- `--sip-degree N` — SIP polynomial degree for gwcs → FITS WCS (default 4, accurate to ~0.1 pixel across a Roman SCA)
- `--no-metadata`, `--metadata-dir` — control the per-exposure metadata CSV
- `--no-catalog`, `--catalog-radius`, `--catalog-color`, `--catalog-extended-color`, `--catalog-label`, `--catalog-show-labels`, `--catalog-include-flagged`, `--catalog-dir`
- `--ds9-target NAME` — DS9 XPA target name
- `--no-dq-overlay` — skip the DQ bad-pixel overlay in DS9
- `--workers N` — concurrent SCA streams (default 8)
- `--coadd` — stream coadd mosaic tiles instead of per-SCA exposures

### Python API

```python
from roman_fits import to_fits_files, to_ds9, download_catalogs
from roman_mast import list_data

res = list_data(program=1039, execution_plan=2, pass_=4, observation=7, sca_only=True)
exposure = res.exposures[0]

# Stream materialized data models {sca: DataModel}
af_dict = res.stream(1)

# Write per-SCA FITS files
to_fits_files(af_dict, exposure, out_dir='output/', compress=False)

# Or pipe every SCA into DS9 as one WCS mosaic
to_ds9(af_dict, exposure, dq_overlay=True)

# Fetch L4 per-SCA catalog parquets for overlay
catalog_paths = download_catalogs(exposure, missions=res.missions)
to_ds9(af_dict, exposure, catalog_paths=catalog_paths)
```

### SIP WCS notes

Roman ships a `gwcs` object per SCA. DS9 / classic astropy consume FITS
headers only, so we approximate `gwcs → FITS SIP` via
`wcs.to_fits_sip(degree=sip_degree)`. The default `degree=4` is accurate
to ~0.1 pixel across a Roman SCA — plenty for DS9 mosaicking and typical
aperture photometry.

---

## `roman-view-sca` — Stream and Visualize a Single SCA

Thin wrapper around `roman_mast` + `roman_fits` for interactive single-SCA
exploration, with optional 2D polynomial background subtraction and
aperture photometry.

### Input modes

1. **Direct S3 or local path** — pass positional `uri` + `filename`.
2. **MAST query** — pass `--program` / `--pass` / etc. plus `--sca N`.

### Display modes

- `--display ds9` — overlay in DS9 (DS9 + XPA must be running). Default.
- `--display mpl` — matplotlib figure (headless-friendly).

### Examples

```bash
# Stream WFI11 from the public tutorial dataset (anonymous S3)
roman-view-sca --display mpl \
    s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
    r0003201001001001004_0001_wfi11_f106_cal.asdf

# Local cache
roman-view-sca --display mpl \
    ./cache r0003201001001001004_0001_wfi11_f106_cal.asdf

# Same, with 2D polynomial background subtraction (model + residuals figure)
roman-view-sca --display mpl --bkg \
    s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
    r0003201001001001004_0001_wfi11_f106_cal.asdf

# Background + aperture photometry together (shares the same polynomial fit)
roman-view-sca --display mpl --bkg --phot \
    s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
    r0003201001001001004_0001_wfi11_f106_cal.asdf

# MAST query mode — pick SCA 4 of exposure 1
roman-view-sca --program 1039 --execution-plan 2 --pass 4 --observation 7 --exposure 1 --sca 4 --display ds9

# Connect to a specific DS9 by XPA target name
roman-view-sca --program 1039 --execution-plan 2 --pass 4 --observation 7 --exposure 1 --sca 1 --ds9 myds9

# Save a matplotlib PNG (headless)
roman-view-sca --display mpl --phot-out phot.csv \
    ./cache r0003201001001001004_0001_wfi11_f106_cal.asdf
```

### Photometry parameters (`--phot`)

| Flag | Default | Meaning |
| --- | --- | --- |
| `--fwhm N` | 1.5 | Source FWHM in pixels |
| `--det-sigma N` | 10.0 | DAOStarFinder detection threshold σ |
| `--aper N` | 2.5 | Aperture radius (× FWHM) |
| `--annulus-inner N` | 6.0 | Inner sky annulus radius (× FWHM) |
| `--annulus-outer N` | 8.0 | Outer sky annulus radius (× FWHM) |
| `--snr-threshold N` | 5.0 | Minimum SNR to keep a source |
| `--bkg-poly-degree N` | 3 | 2D polynomial degree for background fit |
| `--phot-out PATH` | — | Write photometry table to a CSV |

Display polish: `--channels` (H4RG 32-channel grid), `--grid` (8×8 section
grid + labels), `--scale zscale`, `--cmap viridis`, `--sip-degree 4`,
`--no-dq`.

---

## `roman-phot` — Batch Aperture Photometry Across All 18 SCAs

Streams every SCA of an exposure in parallel, runs source detection +
aperture photometry with `photutils`, and writes a combined per-source CSV
plus a per-SCA summary. Optional focal-plane mosaics (background + source
density + image thumbnails).

### Input modes

- **MAST query** — same filters as `roman-mast` (`--program`, `--pass`,
  `--exposures`, `--scas`, …). Forces `--sca-only`.
- **URI file** — plain text, one S3 (or local) URI per line, `#` comments OK.
- **`--remake-mosaics <mosaic_data.npz>`** — regenerate PNG mosaics from a
  previous run without re-doing photometry.

### Examples

```bash
# One line per S3 URI:
#   s3://stpubdata/roman/.../r0103900004001001007_0001_wfi01_f106_cal.asdf
#   s3://stpubdata/roman/.../r0103900004001001007_0001_wfi02_f106_cal.asdf
#   ...
roman-phot --uri-file my_exposure.txt

# Or query MAST directly (no URI file needed)
roman-phot --program 1039 --execution-plan 2 --pass 4 --observation 7 --exposures 1

# Photometry tuning
roman-phot --uri-file my_exposure.txt \
    --fwhm 1.5 --det-sigma 10 --snr-threshold 5

# Also emit per-SCA CSVs
roman-phot --uri-file my_exposure.txt --per-sca

# Background + source-density + image mosaics
roman-phot --uri-file my_exposure.txt --bkg-mosaic --image-mosaic

# Re-render mosaics from a saved .npz (no re-run of photometry)
roman-phot --remake-mosaics r0103900004001001007_0001/mosaic_data.npz
```

### Outputs

All outputs land in `{visit_id}_{exposure_num}/`:

| File | When | Contents |
| --- | --- | --- |
| `sources.csv` | always (if any sources) | Combined per-source table (all SCAs) with `sca`, `detector`, and photutils columns |
| `summary.csv` | always | Per-SCA statistics: `bkg_level`, `bkg_rms`, `n_sources`, … |
| `histograms.png` | unless `--no-hist` | Histograms of `aperture_sum_0` and `snr` |
| `sca{NN}.csv` | with `--per-sca` | Per-SCA source table |
| `bkg_mosaic.png` | with `--bkg-mosaic` | Focal-plane background map |
| `bkg_mosaic_full.png` | with `--bkg-mosaic` | Fine-resolution fitted-background mosaic |
| `source_mosaic.png` | with `--bkg-mosaic` | Focal-plane source density |
| `mosaic_data.npz` | with `--bkg-mosaic` | Raw arrays for `--remake-mosaics` |
| `image_mosaic.png` | with `--image-mosaic` | Focal-plane image thumbnails |
| `sca{NN}.png` | with MAST mode + `--bkg`/`--phot` | Per-SCA image with sources overlaid |

### Photometry parameters

Same defaults as `roman-view-sca` — see the table above.

Background-mosaic tuning: `--bkg-superpixel N`, `--bkg-superpixel-full N`,
`--bkg-mask-sigma N`, `--bkg-dilate N`, `--bkg-poly N`.

Concurrency: `--workers N` (default 8).

---

## `roman-metadata` — Export ASDF Metadata to CSV

Streams the same exposures `roman-fits` operates on, flattens each SCA's
`meta` tree to dot-notation columns, and writes **one CSV row per (exposure,
SCA)**. It does *not* pretty-print a metadata tree — the output is a
spreadsheet-friendly CSV.

`roman-fits` also calls this automatically to drop
`metadata_{visit_id}_exp{NN}.csv` next to each exposure's FITS/DS9 output,
unless `--no-metadata` is passed.

### Examples

```bash
# List matching exposures without writing anything
roman-metadata --program 1039 --execution-plan 2 --pass 4 --observation 7 --list

# Export metadata for exposure 1 (default: --no-list needed to actually write)
roman-metadata --program 1039 --execution-plan 2 --pass 4 --observation 7 --exposures 1 --no-list

# Only some SCAs across a range of exposures
roman-metadata --program 1039 --execution-plan 2 --pass 4 --observation 7 \
    --exposures 1-3 --scas 1-6 --no-list

# Level-1 (uncal) metadata for one visit
roman-metadata --visit-id 0103900004001001007 \
    --data-level 1 --exposures all --no-list

# Custom output path
roman-metadata --program 1039 --execution-plan 2 --pass 4 --observation 7 --exposures 1 \
    --no-list --output /tmp/meta.csv
```

The CSV includes:

- Identifiers: `visit_id`, `exposure`, `sca`
- Data stats: `data_shape`, `data_dtype`, `data_min`, `data_max`,
  `data_mean`, `data_valid_pixels`, `data_nan_pixels`
- Every scalar field under `meta.*`, flattened as `meta.pointing.target_ra`,
  `meta.exposure.effective_exposure_time`, `meta.rcs.led`, …

Rows are sorted by `(visit_id, meta.exposure.start_time, exposure, sca)`
and a warning is emitted if fewer exposures are exported than
`meta.visit.nexposures` claims.

### Python API

```python
from roman_mast import list_data
from roman_metadata import export_csv, extract_rows, write_metadata_csv

res = list_data(program=1039, execution_plan=2, pass_=4, observation=7, sca_only=True)

# One-shot for a range of exposures
export_csv(res, indices=[1, 2, 3], scas=None, output='meta.csv')

# Or reuse an already-streamed exposure (no re-streaming)
af_dict = res.stream(1)
rows = extract_rows(af_dict, res.select(1))
```

---

## `roman-color` — RGB Composite of One SCA in DS9

Builds an RGB composite from three exposures (one per filter) of the same
SCA, aligns them in pixel space, writes the aligned FITS files to disk,
and pushes the trio into DS9 as an RGB frame with auto-computed asinh
limits per channel.

Because Roman commissioning WCS isn't fully calibrated yet, alignment is
done in **raw pixel space** (not WCS reprojection): the three exposures
of the same SCA share a 4088×4088 pixel grid, and any small dither shows
up as an integer-ish (dy, dx) shift measured directly by FFT
cross-correlation.

### Command line

Each layer takes a `key=val,...` spec with any `roman_mast.list_data`
filter (`program`, `pass`, `execution_plan`, `observation`, `visit`,
`optical_element` / `filter` shortcut, `detector`, ...) plus a required
`exposure=N` selecting which exposure inside that query. `--program` and
`--pass` at the top level are inherited by all three layers unless a
layer's own spec overrides them.

```bash
# Select each channel by filter name — program/pass shared globally.
roman-color --sca 3 --program 1047 --pass 1 \
    --blue  filter=F087,exposure=2 \
    --green filter=F129,exposure=2 \
    --red   filter=F184,exposure=2

# Or by observation number, still inheriting the global program/pass.
roman-color --sca 7 --program 1047 --pass 1 \
    --blue  observation=12,exposure=2 \
    --green observation=5,exposure=2 \
    --red   observation=9,exposure=2

# Any layer can override the globals — mix programs/passes across layers.
roman-color --sca 3 --program 1047 --pass 1 \
    --blue  filter=F087,exposure=2 \
    --green filter=F129,exposure=2 \
    --red   program=1052,pass=3,filter=F184,exposure=1

# Reload from the cached aligned FITS in --out-dir without re-streaming
# (great for iterating on stretch / band ratios).
roman-color --sca 3 --program 1047 --pass 1 \
    --blue  filter=F087,exposure=2 \
    --green filter=F129,exposure=2 \
    --red   filter=F184,exposure=2 \
    --from-cache

# Full 18-SCA color mosaic. Streams every SCA in each filter in parallel
# (--workers per filter, default 8), aligns each SCA independently, writes
# 54 aligned FITS to --out-dir, and loads 18 RGB frames into DS9 locked by
# WCS with a shared per-channel asinh stretch.
roman-color --mosaic --program 1047 --pass 1 --workers 8 \
    --blue  filter=F087,exposure=2 \
    --green filter=F129,exposure=2 \
    --red   filter=F184,exposure=2

# --mosaic --from-cache re-uses the 54 aligned FITS without re-streaming.
roman-color --mosaic --program 1047 --pass 1 --from-cache \
    --blue  filter=F087,exposure=2 \
    --green filter=F129,exposure=2 \
    --red   filter=F184,exposure=2
```

### What it does, step by step

1. **Parse the three channel specs.** Each of `--red/--green/--blue` is a
   `key=val,...` string; keys map to `roman_mast.list_data` filters, plus
   a required `exposure=N`.  `--program` and `--pass` at the top level
   inherit into every layer that doesn't already set them.
2. **Query MAST per channel.** For each, calls `list_data(**filters)`
   and picks the exposure whose `exposure` field matches. Warns if
   multiple observations match; add `observation=N` in the spec to
   disambiguate.
3. **Stream one SCA per channel.** Uses
   `roman_fits.stream_materialized(exposure, missions, scas=[sca])` to
   read the ASDF from S3 and materialize the `data` array into memory
   while the pre-signed URL is fresh, then closes the AsdfFile.
4. **Extract data + WCS header.** Pulls `dm.data` (float32) and a SIP
   FITS header from `dm.meta.wcs.to_fits_sip(degree=4)`. WCS goes into
   the FITS output but is **not** used for alignment.
5. **Blue is the alignment pivot.** Blue's raw data is used unchanged as
   the reference; green and red are aligned to it.
6. **Bounded FFT cross-correlation.** For green and red, subtracts the
   sigma-clipped median, zero-fills NaN pixels, computes a full 2D
   `fftconvolve`-based correlation, and picks the peak within a **±10
   px** search window around zero shift. Small window because real
   dither residuals should be sub-pixel to a few px; anything larger is
   spurious lock-on to cross-filter source-density differences.
7. **Apply the (dy, dx) shift.** Zero-fills NaN → `scipy.ndimage.shift`
   with `order=3` → re-NaNs the outside-frame footprint via a
   nearest-neighbour-shifted finite mask (avoids cubic-spline NaN
   smearing across the whole image).
8. **Write per-SCA aligned FITS.** Three files
   (`{channel}_obs{OO}_exp{NN}_{filter}_sca{SS}.fits`) in `--out-dir`
   (default `rgb_sca{NN}/`), all sharing blue's WCS header. These are
   the artifacts you keep for offline band-ratio experiments.
9. **Push to DS9 via pyds9.** `frame new rgb`, then for each channel
   pipes the in-memory FITS bytes to `rgb channel <c>` + `fits`, and
   sets an auto-computed asinh stretch:
   `limits = [median − 1σ, median + gain·σ]` where
   `gain = {blue: 30, green: 25, red: 20}` — F184 gets pushed harder to
   compensate for its lower throughput and keep the composite from
   going blue-dominated.

`--from-cache` skips steps 2–8 entirely: reads the three FITS from disk
(glob-matches by channel + SCA if the exact filename doesn't hit) and
jumps to the DS9 push. The alignment was baked in at cache-write time.

### `--mosaic` mode

`--mosaic` replaces `--sca N` and does the whole focal plane:

1. Calls `stream_materialized(exp, scas=1..18, max_workers=W)` **three
   times** (once per filter). Each call already fans out over `W`
   threads (default 8), so the peak concurrency is 8 streams per filter,
   running sequentially per filter.
2. Loops over SCAs 1..18. For each SCA, runs the same bounded-FFT
   pixel-space alignment used in single-SCA mode (blue as pivot,
   ±10 px search).
3. Writes 54 per-SCA aligned FITS to `--out-dir` (default `rgb_mosaic/`).
4. Computes **one** `(lo, hi)` stretch per channel by pooling
   sigma-clipped stats across all 18 SCAs — so every RGB frame in the
   mosaic gets the same asinh scale, and colors are consistent across
   the focal plane.
5. Builds **three multi-extension FITS** (one per RGB channel, each
   holding all 18 SCAs as extensions with SIP WCS headers) and writes
   them to `out_dir/mosaic_{red,green,blue}.fits`.
6. Pushes the three MEFs into **a single DS9 RGB frame** via
   `fits mosaicimage wcs` on the matching rgb channel — DS9 stitches the
   18 SCAs per channel by WCS and combines the three channels into one
   focal-plane color view.

`--mosaic --from-cache` reloads the 54 FITS from `--out-dir` and jumps
straight to step 5 — fastest way to iterate on the stretch after the
initial run.

### Tuning the stretch

Per-channel σ-gains live at the top of `_run_ds9` in `roman_color.py`:

```python
GAIN_SIGMA = {'blue': 30.0, 'green': 25.0, 'red': 20.0}
LOW_SIGMA = 1.0
```

Bump a channel's gain up if that color is too weak in the composite,
down if it's saturating. `LOW_SIGMA` controls the low cutoff (kills sky
background); raise for a darker sky, lower for more shot noise showing.

---

## `roman-telem` — Query Roman Telemetry from MAST Engineering DB

Query telemetry mnemonics with optional CSV/Parquet/HDF5/Pickle output
and optional trending plots via `roman_telem_plot`.

### Command line

```bash
# Query a few mnemonics directly
roman-telem --mnemonics WFI_MCE_SRCS_PD1_V WFI_MCE_SRCS_PD2_V \
    -s "2026-05-01" -e "2026-07-02" \
    --output data.csv

# Query mnemonics listed in a text file (whitespace, newlines, commas OK; # comments OK)
roman-telem --mnemonics-file mnemonics.txt \
    -s "2025-01-01" -e "2025-01-02" \
    --output data.parquet

# Discover available groups in a YAML plot-groups config
roman-telem --plot-groups tlm_groups.yaml --list-groups

# Query a subset of groups (much faster than "all")
roman-telem --plot-groups tlm_groups.yaml \
    --select-groups "rcs_pd" \
    -s "2026-05-01" -e "2026-07-02"

# Query, save data, AND produce a trending plot in one shot
roman-telem --plot-groups tlm_groups.yaml \
    --select-groups "rcs_pd" \
    -s "2026-05-01" -e "2026-07-02" \
    --output rcs_pd.csv \
    --plot-output rcs_pd_trend.png \
    --plot-layout vertical
```

### Options

**Mnemonic input** (pick any combination — results are de-duplicated):
- `--mnemonics NAME [NAME ...]` — one or more mnemonics
- `--mnemonics-file PATH` — text file, one mnemonic per line (whitespace / comma separated, `#` comments OK)
- `--plot-groups FILE.yaml` — load mnemonics from a plot-groups YAML

**Time range** (required unless `--list-groups`):
- `-s`, `--start` — ISO datetime, e.g. `"2026-05-01T00:00:00"`
- `-e`, `--end`

**Filtering / discovery**:
- `--select-groups "g1,g2"` — restrict to these groups from `--plot-groups`
- `--list-groups` — print groups (name / label / mnemonic count) and exit; requires `--plot-groups`

**Output**:
- `--output FILE` — save to CSV / Parquet / HDF5 / Pickle (format inferred from extension)
- `--output-format {csv,parquet,hdf5,pickle}` — override the inferred format
- `--no-combine` — write one file per mnemonic (`base_MNEMONIC.ext`) instead of a combined DataFrame

**Server**:
- `--api-token TOKEN` — override the env-var token

**Retry / verbosity**:
- `--max-retries N` (default 3), `--retry-delay S` (default 1.0)
- `--quiet`

**Plotting** (also triggers a matplotlib plot):
- `--plot-output PATH` — save the trending plot to PNG/PDF
- `--plot-layout {vertical,horizontal,grid}` (default `vertical`)
- `--show` — display interactively

### Plot groups

Plot groups are defined in a YAML config. Example (from the shipped
`tlm_groups.yaml`):

```yaml
groups:
  rcs_pd:
    label: "RCS photodiodes"
    y_label: "Voltage (V)"
    mnemonics:
      - WFI_MCE_SRCS_PD1_V
      - WFI_MCE_SRCS_PD2_V
```

Structure:

```yaml
groups:
  <group_name>:
    label: <string, shown on subplot title>
    y_label: <string, shown on y-axis>       # optional
    mnemonics:
      - <MNEMONIC_NAME>
      - <MNEMONIC_NAME>
      ...
```

### Python API

```python
from roman_telem import MASTEngDBQuery, query_telemetry

# Low-level: use the query class directly
query = MASTEngDBQuery(api_token="your_token", server="mast")

df = query.query_mnemonic(
    "WFI_MCE_SRCS_PD1_V",
    start_time="2026-05-01",
    end_time="2026-05-02",
)

dfs = query.query_multiple_mnemonics(
    ["WFI_MCE_SRCS_PD1_V", "WFI_MCE_SRCS_PD2_V"],
    start_time="2026-05-01",
    end_time="2026-05-02",
    combine=True,  # → single wide DataFrame indexed by ObsTime
)

# High-level: same convenience the CLI uses
df = query_telemetry(
    groups_config="tlm_groups.yaml",
    selected_groups=["rcs_pd"],
    start_time="2026-05-01",
    end_time="2026-05-02",

    combine=True,
)
```

---

## `roman-telem-plot` — Plot Telemetry Data

Standalone plotting for telemetry data written by `roman-telem`. Same
underlying `plot_telemetry()` function that `roman-telem --plot-output`
uses.

### Command line

```bash
# Plot every mnemonic present in the file
python -m roman_telem_plot data.csv --show

# Only selected groups, with a suptitle
python -m roman_telem_plot data.csv \
    --plot-groups tlm_groups.yaml \
    --select-groups "rcs_pd" \
    --suptitle "RCS photodiodes trending" \
    --layout vertical \
    --output trending.png --dpi 300
```

> The tool is exposed as `roman-telem-plot` if the entry point is
> installed; equivalently `python -m roman_telem_plot` always works.

### Options

**Input** — a positional `INPUT` file: CSV, Parquet, HDF5, or Pickle
produced by `roman-telem`.

**Grouping / filtering**:
- `--plot-groups FILE.yaml` — load groups from a YAML file (subplot per group)
- `--select-groups "g1,g2"` — restrict to these groups
- `--mnemonics "M1,M2"` — restrict to these mnemonics (no grouping)

**Display**:
- `--layout {vertical,horizontal,grid}` (default `vertical`)
- `--suptitle "Title"` — figure-level title
- `--output FILE` — save to PNG/PDF (default: show interactively)
- `--dpi N` (default 120)
- `--time-col NAME`, `--value-col NAME` — override auto-detected columns

### Python API

```python
import pandas as pd
from roman_telem_plot import plot_telemetry, plot_grouped, plot_mnemonics

df = pd.read_csv("data.csv")   # or read_parquet / read_hdf / read_pickle

# Grouped subplots
plot_telemetry(
    df,
    groups_config="tlm_groups.yaml",
    selected_groups=["rcs_pd"],
    layout="horizontal",
    output="trending.png",
    show=False,
)

# All mnemonics on a single Axes
plot_mnemonics(df, output="all.png")
```

The input DataFrame is expected in **long format** (mnemonic column +
time column + value column). Column names are auto-detected from common
choices (`mnemonic`, `ObsTime` / `time` / `datetime`, `EUValue` / `value`
/ `eng_value`, …) — override with `--time-col` / `--value-col` if needed.

---

## CAR (Commissioning Activity Report) Overlays on Telemetry Plots

Telemetry plots can be annotated with **CAR activity spans** — colored vertical shaded regions labeled with CAR numbers and activity names. This helps correlate operations (thruster burns, slews, deployments, etc.) with telemetry anomalies or trending changes.

CAR data comes from the **Roman commissioning operations CAST (Commissioning Activity Spreadsheet Template)**, which is publicly archived at STScI.

### Workflow: Download CAST and Extract CAR Summary

#### Step 1: Download the CAST from STScI

1. Visit the **[Roman Mission Operations Archive](https://www.stsci.edu/roman/mission-operations/)** (or search for "Roman CAST as-run")
2. Download the **RST_Commissioning_CAST_AS_RUN.xlsx** file (or current version; filename may vary by revision)
3. Save it locally, e.g. to your working directory

#### Step 2: Extract CAR information using `extract_car_info.py`

The `extract_car_info.py` script parses the CAST Excel file and extracts commissioning activities:

```bash
python extract_car_info.py RST_Commissioning_CAST_AS_RUN.xlsx
```

**Output:**
- A `RST_Commissioning_CAST_AS_RUN_CAR_Summary.csv` file with columns:
  - `CAR_Number` (e.g., `CAR-173.5`)
  - `CAR_Name` (e.g., `Count-rate Dependent Nonlinearity (Direct)`)
  - `APT_Program` (e.g., `1025`)
  - `MET` (Mission Elapsed Time, e.g., `L + 018:11:12:00`)
  - `Start_Time` (wall-clock time, UTC)
  - `Duration` (e.g., `03:10:00`)

**Example:**
```
CAR_Number,CAR_Name,APT_Program,MET,Start_Time,Duration
CAR-173.5,Count-rate Dependent Nonlinearity (Direct),1025,L + 018:11:12:00,19:56:00,03:10:00
CAR-190.20,RCS Daily Monitor,1039,L + 016:21:18:00,09:18:00,00:32:00
CAR-086.4,ST-FGS Alignment - Confirmation Image & WIM Closed-loop Checkout,1028,L + 017:19:56:00,03:06:00,02:10:00
...
```

**Features:**
- Automatically skips CARs marked as "Crossed Out" (strikethrough formatting in the Excel file)
- Extracts APT program IDs from both CAR headers and step descriptions
- Handles malformed or missing time fields gracefully

#### Step 3: Place `car_summary.csv` in the tools directory

Rename or symlink the extracted CSV to `car_summary.csv` in the `roman-mast-tools` directory:

```bash
cp RST_Commissioning_CAST_AS_RUN_CAR_Summary.csv roman-mast-tools/car_summary.csv
```

Once present, `roman-telem` will automatically detect and use it for plot annotations.

### Using CAR Overlays with `roman-telem`

#### Basic usage: Plot telemetry with CAR activity bands

```bash
# Query temperature data for a subsystem during a CAR-heavy period
roman-telem --plot-groups tlm_groups.yaml \
    --select-groups "temp_fps" \
    -s "2026-09-16" -e "2026-09-20" \
    --show
```

**What you'll see:**
1. **Three subplots** (one per temperature sensor group in `temp_fps`):
   - Each subplot shows temperature trending over time
   - Colored vertical shaded bands indicate CAR activity windows
   - Color is consistent within a program (e.g., all Program 1025 activities are blue)

2. **Leader-line callouts** below the bottom axis:
   - Vertical stubs drop from each CAR's midpoint
   - Text labels (e.g., `CAR-173.5/1025: Count-rate Dependent Nonlinearity…`) tilt at 45° pointing upward into the stub
   - Multiple tiers of labels automatically deconflict overlapping CARs
   - Labels are truncated to 50 characters with an ellipsis

3. **Date/time labels** on the x-axis:
   - Format: `259|09-16 12` (day-of-year | month-day hour)
   - Tilted 45° to match the visual language of the callouts

#### Example: Correlate FPS temperatures with CAR activities

During commission, FPS (Focal Plane Subsystem) temperatures fluctuate due to operational activities. The CAR overlay helps identify which activity caused each spike:

```bash
# Plot FPS temperature groups for the entire first week of operation
roman-telem --plot-groups tlm_groups.yaml \
    --select-groups "temp_fps" \
    -s "2026-09-10" -e "2026-09-17" \
    --output fps_temps_with_cars.png \
    --plot-layout vertical
```

**What this reveals:**
- `CAR-078.1: PAM ReCalibration` (P1104) on 2026-09-10 shows a sustained temperature rise
- `CAR-190.x: RCS Daily Monitor` (P1039) series on 2026-09-12–09-16 produce brief spikes
- Temperature returns to nominal after `CAR-086.4: ST-FGS Alignment` completes on 2026-09-18

This makes it easy to debug thermal issues or validate that temperatures are responding as expected to each activity.

#### Command-line options for CAR overlays

- `--program-spans` — (legacy) query MAST for Roman exposure windows instead of loading CAR CSV. CARs take precedence if `car_summary.csv` is found.
- `--car-csv PATH` — specify a custom CAR CSV path (default: looks for `car_summary.csv` in the package directory)

### Python API: Plot with CAR Overlays

```python
import pandas as pd
from roman_telem_plot import plot_telemetry
from roman_telem_cars import load_car_spans, ROMAN_LAUNCH

# Load telemetry data
df = pd.read_csv("temp_fps.csv")

# Load CAR spans from the extracted CSV
car_spans = load_car_spans(
    "car_summary.csv",
    launch_time=ROMAN_LAUNCH,
    start_window=pd.Timestamp("2026-09-16"),
    end_window=pd.Timestamp("2026-09-20"),
)

# Plot with CAR overlay
plot_telemetry(
    df,
    groups_config="tlm_groups.yaml",
    selected_groups=["temp_fps"],
    layout="vertical",
    program_spans=car_spans,
    output="fps_temps_with_cars.png",
    show=False,
)
```

### Interpreting the CAR Label Format

Each callout label follows the pattern:

```
CAR-173.5/1025: Count-rate Dependent Nonlinearity…
├──────┬─────┼────┬──────────────────────────────┘
  CAR    subID APT  Activity name (truncated at 50 chars)
```

- **CAR number** (e.g., `CAR-173.5`): unique activity identifier from the CAST
- **APT program** (e.g., `1025`): science program that scheduled this activity
- **Activity name**: human-readable description (used for visual correlation with telemetry trends)

The label color matches the program's shaded bands — follow color → shaded region → telemetry effect visually.

### Troubleshooting

| Issue | Solution |
|-------|----------|
| CAR bands not showing on plot | Verify `car_summary.csv` exists in the package directory; check that date range overlaps CAR times |
| Labels too crowded or missing | Increase `--days` range or reduce subplot count; very dense days may exceed 20 tiers (hard safety cap) |
| Callouts overlapping with date labels | This is rare with constrained layout; if it happens, the auto-sizing should resolve it on next run |
| MET times don't match wall-clock times | Check ROMAN_LAUNCH constant in `roman_telem_cars.py`; it should match the mission epoch used in your CAST |

---

## Architecture

### Module hierarchy

**Science data pipeline**
```
roman_mast.py       ← foundation: MAST auth, product search, S3 streaming
roman_fits.py       ← output layer: ASDF → FITS/DS9, SIP WCS, catalog overlays
roman_metadata.py   ← metadata CSV export
roman_view_sca.py   ← interactive single-SCA viewer + background/photometry
roman_phot.py       ← batch photometry across all 18 SCAs + focal-plane mosaics
```

**Telemetry pipeline**
```
roman_telem.py      ← MAST EDB queries: mnemonic search, retry logic, I/O
roman_telem_plot.py ← trending plots with grouping + layout control
```

Every science tool builds on `roman_mast` primitives. Every telemetry tool
builds on `roman_telem`'s query and I/O logic.

### Streaming pipeline

1. **Query** — MAST search returns a `DataResults` with grouped `Exposure`
   objects (per `(visit_id, exposure)`).
2. **Stream** — anonymous S3 access via s3fs, or authenticated MAST API
   via astroquery's `MastMissions.read_product`.
3. **Load** — `roman_datamodels` + ASDF parsing into memory.
4. **WCS** — SIP approximation computed by `roman_fits` for FITS
   compatibility (`gwcs → wcs.to_fits_sip(degree=sip_degree)`).
5. **Dispatch** — FITS files, DS9 mosaic, Imviz, in-memory analysis, or
   metadata CSV.

### Focal-plane layout

The 18 WFI SCAs are arranged in a 3×6 grid. Rotations (per
`roman_view_sca._SCA_ROTATION`):

- **No rotation:** SCAs 3, 6, 9, 12, 15, 18 (every third SCA — the `r=0` row)
- **180° rotation:** all others

Background mosaics, DS9 displays, and matplotlib plots respect this layout.

---

## Dependencies

**Core**
- `numpy`, `astropy` — array math, WCS, FITS I/O, tables
- `roman_datamodels`, `romancal`, `rad` — Roman-specific ASDF schemas and datamodels
- `asdf` — ASDF file parsing
- `gwcs` — generalized WCS (Roman ships gwcs per SCA; approximated to SIP)
- `photutils` — source detection, aperture photometry, background estimation
- `matplotlib` — visualization
- `s3fs`, `fsspec` — anonymous S3 streaming
- `astroquery` — MAST authenticated queries (`MastMissions`)
- `pandas`, `pyarrow` — telemetry DataFrames, Parquet I/O
- `pyyaml` — plot-group config parsing
- `requests` — MAST Engineering DB REST client
- `keyring`, `python-dotenv` — credential handling

**Optional**
- **DS9 + XPA** — required for `--display ds9` and `roman-fits --ds9`
  (local installation, DS9 must be running before the tool is invoked)
- `h5py` / `tables` — HDF5 output for `roman-telem`
- `tqdm` — progress bars where available

See [Installation](#installation) above.

---

## Configuration

### Output directory

`roman-phot` (and any tool that writes per-exposure folders) reads an
optional `config.ini` to set the output root:

```ini
# ~/config.ini  or  ./config.ini
[paths]
output_root = ~/roman_output
```

Search order:

1. Path given via `--config` (if the tool supports it)
2. `~/config.ini`
3. `./config.ini`
4. Fallback: current working directory

Environment variables (`~`, `$HOME`, etc.) are expanded.

### Environment variables

| Variable | Purpose |
|---|---|
| `MAST_API_TOKEN` | MAST auth token (`mast.stsci.edu`) |
| `AWS_PROFILE` | Optional; anonymous S3 works without any AWS credentials |

Tokens can also be passed explicitly via `--api-token` on the CLI or the
`api_token=` kwarg in the Python API.

---

## Tips and Gotchas

**Streaming vs. caching.** Streaming from S3 is convenient but slow for
iterative work. For repeated analysis, download once with
`bash download_all_scas.sh` and pass local paths.

**Concurrency.** `roman-fits`, `roman-metadata`, and `roman-phot` all
support `--workers N` for concurrent SCA streams. Default is 8; set to 1
for sequential mode when debugging.

**Memory.** A full exposure (18 SCAs, level-2 `cal`) is roughly 1.2 GB of
`data` + `dq` in memory. `roman_fits` materializes `data` and `dq` only;
add `err` / `var_*` to `_MATERIALIZE_ATTRS` in `roman_fits.py` if you
need them (each adds ~17-67 MB per SCA × 18 SCAs).

**SIP degree.** The default `sip_degree=4` is accurate to ~0.1 pixel
across a Roman SCA. Bump to 5 for sub-pixel precision at the corners,
or drop to 3 for faster export if you don't care about mosaic edges.

**DS9 not showing up.** DS9 must already be running before you invoke
`--display ds9` or `roman-fits --ds9`. If you have multiple DS9 windows
open, pick one with `--ds9-target <name>`.

**Catalog overlays.** Not every SCA has an L4 catalog on MAST — those
SCAs silently overlay nothing. `roman-fits` writes one parquet per SCA
that has a catalog plus one combined `.reg` file DS9 loads.

**Telemetry time windows.** MAST EDB queries are inclusive on both ends
and time strings are interpreted as UTC. Very long windows (weeks to
months) on high-rate mnemonics can return millions of samples — use
`--select-groups` and short windows when iterating.

**"No mnemonics specified".** `roman-telem` needs at least one of
`--mnemonics`, `--mnemonics-file`, or `--plot-groups`. Combining them is
fine — the union is de-duplicated.

---

## FAQ

**Q: Can I use this without MAST authentication?**
A: Yes, for all public tutorial data via anonymous S3
(`s3://stpubdata/...`). Authentication is only needed for private
programs and any telemetry (`roman-telem`, `roman-telem-plot`).

**Q: How do I cache files locally for fast iteration?**
A: Run `bash download_all_scas.sh` to pull all 18 tutorial SCAs (~3.4 GB)
into `./cache/`, then use `--uri ./cache --filename ...` or pass local
paths in a URI file.

**Q: Can I run this headless / in a container?**
A: Yes. Use `--display mpl` for matplotlib output. DS9 display requires
a running DS9 instance and XPA; use X11 forwarding or a virtual display
(`xvfb-run`) if you need DS9 in a container.

**Q: What does the "source mask" percentage mean in background maps?**
A: The fraction of pixels masked during background fitting (bright
sources, bad pixels, cosmic rays, etc.). Values of 15–20% are typical
for science data.

**Q: Do I need a MAST token for telemetry queries?**
A: Yes. The MAST Engineering Database requires authentication. Obtain a
token from the MAST Portal (Settings → API Tokens) and store it in
`MAST_API_TOKEN`.

**Q: How do I use plot groups for telemetry?**
A: Define groups in a YAML config file (e.g. `tlm_groups.yaml`), then
pass it to `roman-telem` with `--plot-groups`. Use `--list-groups` to
discover available names and `--select-groups "g1,g2"` to filter.

**Q: My query is slow / times out — what can I do?**
A: Narrow the time window, restrict to a specific group with
`--select-groups`, and bump `--max-retries` / `--retry-delay`. Save
results to Parquet (`--output data.parquet`) rather than CSV for large
returns.

**Q: The SIP WCS looks off at the edges — why?**
A: Roman's true WCS is a gwcs distortion model. `roman_fits` approximates
it as a SIP polynomial (degree 4 by default), which is accurate to
~0.1 pixel across the SCA. For sub-pixel edge precision, increase the
degree; for cross-SCA astrometry, use the native gwcs directly from the
ASDF file.

---

## References

- [Roman Space Telescope Documentation](https://roman.gsfc.nasa.gov/)
- [Roman User Documentation (RDox)](https://roman-docs.stsci.edu/)
- [MAST Archive](https://mast.stsci.edu/)
- [MAST Engineering Database](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)
- [`roman_datamodels`](https://github.com/spacetelescope/roman_datamodels)
- [`romancal`](https://github.com/spacetelescope/romancal)
- [`photutils` Documentation](https://photutils.readthedocs.io/)
- [`astropy` Documentation](https://docs.astropy.org/)
- [DS9](https://sites.google.com/cfa.harvard.edu/saoimageds9) / [XPA](https://github.com/ericmandel/xpa)

---

## License

See `LICENSE` in the repository root.

## Support

Issues and pull requests welcome. For questions specific to Roman data
products or MAST access, please refer to the Roman Help Desk at STScI.
