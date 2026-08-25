# roman-mast-tools

Tools for querying, streaming, and analyzing Roman WFI data from MAST.

## Installation

Create and activate the conda environment before using any of the tools.
The environment installs all dependencies (including an unreleased `astroquery`
build that adds Roman MAST streaming support) and the package itself.

```bash
conda env create -f environment.yml
conda activate roman-mast-tools
```

If the environment already exists and you need to update it:

```bash
conda env update -f environment.yml --prune
```

All commands below assume the `roman-mast-tools` environment is active.

This installs the following command-line tools:
- `roman-mast` — query and list Roman observations from MAST
- `roman-fits` — convert ASDF to FITS and display in DS9
- `roman-metadata` — extract and display metadata from ASDF files
- `roman-view-sca` — stream and visualize a single SCA with photometry
- `roman-phot` — batch photometry across all SCAs of an exposure

---

## `roman-mast` — Query and Stream from MAST

`roman_mast` is the foundation of the stack. Every other tool builds on its
streaming primitives. It handles MAST authentication, product discovery, and
streaming Roman ASDF files directly into memory.

### Basic Usage

```bash
# List all exposures for a program/pass
roman-mast --program 114 --pass 57 --list

# List only the SCAs (L2 cal files) for a specific detector
roman-mast --program 114 --pass 57 --detector wfi04 --list
```

### Python API

```python
from roman_mast import list_data, print_summary

# Search MAST — all filters are optional
res = list_data(program=114, pass_=57, sca_only=True)
print_summary(res)

# Stream one exposure into memory: {sca: AsdfFile}
af_dict = res.stream('r0011401057001001001_0001')

# One-liner: stream → FITS files on disk
res.to_fits('r0011401057001001001_0001', out_dir='output/')

# One-liner: stream → DS9 mosaic (DS9 must be running)
res.to_ds9('r0011401057001001001_0001')
```

### Filter Reference

| Filter | Type | Description |
|---|---|---|
| `program` | int | APT program ID (e.g. 114) |
| `pass_` | int | Pass within the execution plan |
| `execution_plan` | int | Execution plan within the program |
| `segment` | int | Segment within the pass |
| `observation` | int | Observation within the segment |
| `visit` | int | Visit within the observation |
| `detector` | str/int | `'WFI04'` / `'wfi04'` / `4` |
| `visit_id` | str | Full ID or wildcard, e.g. `'0011401057*'` |
| `exposure` | int | Matches last 4 digits of observation_id |
| `optical_element` | str | e.g. `'F062'`, `'F129'` |
| `exposure_type` | str | e.g. `'WFI_IMAGE'`, `'WFI_DARK'` |
| `product_type` | str | `'l2'` (per-SCA) or `'p_visit_coadd'` (mosaic tiles) |
| `sca_only` | bool | Shortcut for `product_type='l2'`; drops mosaic tiles |
| `data_level` | 1/2/`'gw'`/None | 1→`_uncal`, 2→`_cal`, `'gw'`→`_gw`, None→all |

### Authentication

MAST authentication is only required for private programs. All public data
(including the tutorial dataset) can be accessed anonymously via S3.

```bash
# Store your MAST API token (one-time setup)
export MAST_API_TOKEN=your_token_here
```

---

## `roman-fits` — Convert ASDF to FITS and Display in DS9

`roman_fits` builds on `roman_mast` to convert Roman ASDF files to FITS,
approximating the gwcs WCS as a SIP polynomial for DS9 and standard FITS
tooling compatibility.

### Command Line

```bash
# Convert a local ASDF file to FITS
roman-fits r0003201001001001004_0001_wfi11_f106_cal.asdf

# Convert and open all SCAs in DS9 as a mosaic
roman-fits --ds9 r0003201001001001004_0001_wfi11_f106_cal.asdf
```

Output: `r0003201001001001004_0001_wfi11_f106_cal.fits`

### Python API

```python
from roman_fits import to_fits_files, to_ds9, download_catalogs

# Write one FITS file per SCA
to_fits_files(af_dict, exposure, out_dir='output/', compress=False)

# Pipe all SCAs into a DS9 mosaic via XPA
# Optionally overlay DQ masks and/or L4 catalog sources
to_ds9(af_dict, exposure, dq_overlay=True, catalog_paths=catalog_dict)

# Download per-SCA L4 catalogs (parquet) from MAST
catalog_paths = download_catalogs(exposure, missions=['roman'])
```

### SIP WCS Notes

Roman uses a `gwcs` object per SCA. `roman_fits` approximates this as a
FITS SIP polynomial (`degree=4` by default), accurate to ~0.1 pixel across
a Roman SCA — sufficient for DS9 mosaicking and standard aperture photometry.

---

## `roman-view-sca` — Stream and Visualize a Single SCA

Thin wrapper around `roman_mast` and `roman_fits` for interactive,
single-SCA exploration with optional background subtraction and photometry.

### Input Modes

1. **S3 streaming** (anonymous): pass an S3 base URI and filename
2. **Local files**: pass a local directory and filename
3. **MAST query** (token required): `--program`, `--pass`, `--sca`

### Display Modes

- `--display ds9` — overlay in DS9 (DS9 must be running)
- `--display mpl` — matplotlib figure (headless-friendly)

### Examples

```bash
# Stream WFI11 from the public tutorial dataset
roman-view-sca --display mpl \
  s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
  r0003201001001001004_0001_wfi11_f106_cal.asdf

# Local cache with background subtraction + photometry
roman-view-sca --display mpl --bkg --phot --phot-out wfi11_phot.csv \
  ./cache r0003201001001001004_0001_wfi11_f106_cal.asdf

# MAST query: view SCA 11 from the 2nd exposure with photometry
roman-view-sca --program 114 --pass 57 --exposure 2 --sca 11 \
  --phot --phot-out phot.csv
```

### Analysis Options

- `--bkg` — fit 2D polynomial background; show model + residuals
- `--phot` — run aperture photometry; overlay sources on image
- `--bkg --phot` — combined (shares the same polynomial fit)

Photometry output columns: `id, x_center, y_center, aperture_sum_0, aperture_sum_1, local_bkg_per_pix, local_bkg_total, flux_bkgsub, flux_err, snr`

---

## `roman-phot` — Batch Photometry Across All SCAs

Runs source detection and aperture photometry across all 18 SCAs of an
exposure, assembling results into a combined source catalog and summary
statistics.

### Usage

```bash
# Create a URI file (one S3 path or local path per line)
roman-phot --uri-file my_exposure.txt --bkg-mosaic

# High-sensitivity faint source detection
roman-phot --uri-file my_exposure.txt \
  --fwhm 2.0 --det-sigma 3 --snr-threshold 2 \
  --bkg-mosaic --per-sca
```

### Options

**Background modeling:**
- `--bkg-mosaic` — fit per-SCA background; visualize focal-plane mosaic
- `--bkg-poly N` — polynomial degree (default: 3)
- `--bkg-superpixel N` — superpixel size for background map (default: 8)

**Photometry parameters:**
- `--fwhm N` — expected PSF FWHM in pixels (default: 1.5)
- `--det-sigma N` — detection threshold in sigma (default: 10.0)
- `--snr-threshold N` — keep sources with SNR ≥ N (default: 5.0)
- `--aper N` — aperture radius as multiple of FWHM (default: 2.5)
- `--annulus-inner N`, `--annulus-outer N` — background annulus radii (default: 6.0, 8.0)

### Outputs

- `roman_phot_combined.csv` — source catalog with photometry (all SCAs)
- `roman_phot_summary.csv` — per-SCA statistics
- `roman_phot_bkg_mosaic.png` — background map in WFI focal-plane layout
- `roman_phot_histograms.png` — SNR, background level, and source count distributions
- Per-SCA CSVs (if `--per-sca`)

---

## `roman-metadata` — Inspect ASDF Metadata

```bash
roman-metadata r0003201001001001004_0001_wfi11_f106_cal.asdf
```

Prints the full metadata tree from an ASDF file: observation context, filter,
exposure time, detector, WCS reference, etc.

---

## Architecture

### Module Hierarchy

```
roman_mast.py       ← foundation: MAST auth, product search, S3 streaming
roman_fits.py       ← output layer: ASDF → FITS/DS9, SIP WCS, catalog overlays
photometry.py       ← SourcePhotometry class: detection, aperture phot, PSF phot
roman_view_sca.py   ← interactive single-SCA viewer
roman_phot.py       ← batch photometry across all 18 SCAs
```

Every tool builds on `roman_mast` primitives — streaming is never
re-implemented elsewhere.

### Streaming Pipeline

1. **Query** — MAST search returns a `DataResults` with grouped `Exposure` objects
2. **Stream** — anonymous S3 access via s3fs, or authenticated MAST API
3. **Load** — `roman_datamodels` + ASDF parsing into memory
4. **WCS** — SIP approximation computed by `roman_fits` for FITS compatibility
5. **Dispatch** — FITS files, DS9 mosaic, or in-memory analysis

### Focal-Plane Layout

The 18 WFI SCAs are arranged in a 3×6 grid with specific rotations:
- SCAs 3, 6, 9, 12, 15, 18 — no rotation
- All others — 180° rotation

Background mosaics and DS9 displays respect this layout.

---

## Quick Start: Tutorial Data

### Cache Locally (~3.4 GB)

```bash
bash download_all_scas.sh
```

Creates a `cache/` directory with all 18 tutorial SCAs for fast offline analysis.

### End-to-End: Stream → Photometry → Inspect

```bash
# 1. List available exposures
roman-mast --program 114 --pass 57 --list

# 2. View one SCA interactively
roman-view-sca --display mpl --bkg --phot \
  s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
  r0003201001001001004_0001_wfi11_f106_cal.asdf

# 3. Run batch photometry on all 18 SCAs
roman-phot --uri-file my_exposure.txt --bkg-mosaic
```

---

## Dependencies

Core:
- `numpy`, `astropy` — data structures, WCS
- `roman_datamodels`, `romancal`, `rad` — Roman-specific ASDF formats
- `photutils` — source detection, aperture photometry
- `matplotlib` — visualization
- `s3fs`, `fsspec` — S3 streaming
- `keyring`, `python-dotenv` — credential handling

Optional (for DS9):
- DS9 + XPA (local installation required)

---

## FAQ

**Q: Can I use this without MAST authentication?**

A: Yes. All tutorial examples use anonymous S3 access (`s3://stpubdata/...`). Authentication is only needed for private programs.

**Q: How do I cache files locally for fast iteration?**

A: Run `bash download_all_scas.sh` to download all 18 tutorial SCAs. Then use `--uri ./cache --filename ...` or pass local paths in a URI file.

**Q: Can I run this headless / in a container?**

A: Yes. Use `--display mpl` for matplotlib output. DS9 display requires a running DS9 instance and XPA; use X11 forwarding if needed.

**Q: What does the "source mask" percentage mean in background maps?**

A: The fraction of pixels masked during background fitting (bright sources, bad pixels, etc.). Values of 15–20% are typical for science data.

---

## References

- [Roman Space Telescope Documentation](https://roman.gsfc.nasa.gov)
- [MAST Archive](https://mast.stsci.edu)
- [photutils Documentation](https://photutils.readthedocs.io)
- [astropy Documentation](https://docs.astropy.org)
