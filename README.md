# roman-mast-tools

Tools for streaming Roman WFI data from MAST and analyzing with aperture photometry.

## Installation

```bash
pip install -e .
```

This installs the following command-line tools:
- `roman-mast` — query and list Roman observations from MAST
- `roman-fits` — convert ASDF to FITS
- `roman-metadata` — extract and display metadata
- `roman-view-sca` — stream and visualize a single SCA with photometry
- `roman-phot` — batch photometry across all SCAs of an exposure

## Quick Start

### View a Single SCA

Stream WFI11 from the public tutorial dataset with matplotlib display:

```bash
python roman_view_sca.py --display mpl \
  s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
  r0003201001001001004_0001_wfi11_f106_cal.asdf
```

### Run Photometry on All SCAs

Process all 18 SCAs of an exposure with source detection and aperture photometry:

```bash
# Create a URI file (one S3 path per line)
cat > my_exposure.txt <<EOF
s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/r0003201001001001004_0001_wfi01_f106_cal.asdf
s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/r0003201001001001004_0001_wfi02_f106_cal.asdf
...
s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/r0003201001001001004_0001_wfi18_f106_cal.asdf
EOF

# Run photometry with background subtraction
python roman_phot.py --uri-file my_exposure.txt --bkg-mosaic
```

This generates:
- `roman_phot_combined.csv` — source catalog with photometry (all SCAs)
- `roman_phot_summary.csv` — per-SCA statistics
- `roman_phot_bkg_mosaic.png` — background map in WFI focal-plane layout
- `roman_phot_histograms.png` — SNR and other distributions

## Example Workflows

### 1. Cache Tutorial Data Locally

Download all 18 tutorial SCAs for fast, offline analysis:

```bash
bash download_all_scas.sh
```

This creates a `cache/` directory with all files (~3.4 GB).

### 2. Visualize SCA with Background Subtraction and Photometry

Stream from S3, fit a 2D polynomial background, detect sources, and run aperture photometry:

```bash
python roman_view_sca.py --display mpl --bkg --phot --phot-out wfi11_phot.csv \
  s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
  r0003201001001001004_0001_wfi11_f106_cal.asdf
```

Output: `wfi11_phot.csv` with columns: `id, x_center, y_center, aperture_sum_0, aperture_sum_1, local_bkg_per_pix, local_bkg_total, flux_bkgsub, flux_err, snr`

### 3. Full Exposure Analysis with Custom Photometry Parameters

Detect faint sources with SNR ≥ 3, use 2.0 pixel FWHM:

```bash
python roman_phot.py --uri-file my_exposure.txt \
  --fwhm 2.0 --det-sigma 5 --snr-threshold 3 \
  --bkg-mosaic --per-sca
```

Generates per-SCA CSV files in addition to combined results.

### 4. Batch Process Many Exposures

Script to process a list of exposure URIs:

```bash
for exposure_uri in $(cat exposure_list.txt); do
  echo "Processing: $exposure_uri"
  
  # Create URI file for all 18 SCAs
  cat $(for i in {01..18}; do
    echo "$exposure_uri/r*_wfi${i}_*.asdf" | xargs ls -1 2>/dev/null
  done) > temp_uris.txt
  
  # Run photometry
  python roman_phot.py --uri-file temp_uris.txt --bkg-mosaic
  
  # Archive results
  mkdir -p results/$(basename $exposure_uri)
  mv roman_phot_*.csv roman_phot_*.png results/$(basename $exposure_uri)/
done
```

### 5. Query MAST for Observations (Requires MAST Token)

List all exposures from a program/pass combination:

```bash
python roman_view_sca.py --program 114 --pass 57 --list
```

View SCA 1 from the first matching exposure in DS9:

```bash
python roman_view_sca.py --program 114 --pass 57 --sca 1
```

View SCA 11 from the 2nd exposure with photometry:

```bash
python roman_view_sca.py --program 114 --pass 57 --exposure 2 --sca 11 \
  --phot --phot-out phot.csv
```

### 6. Convert ASDF to FITS

```bash
python roman_fits.py r0003201001001001004_0001_wfi11_f106_cal.asdf
```

Outputs `r0003201001001001004_0001_wfi11_f106_cal.fits`.

## Tools Reference

### `roman_view_sca.py`

Stream and visualize a single SCA with optional photometry overlay.

**Input modes:**
1. **S3 streaming** (anonymous): `--uri s3://... --filename r0003...asdf`
2. **Local files**: `--uri ./cache --filename r0003...asdf`
3. **MAST query** (token required): `--program 114 --pass 57 --sca 1`

**Output modes:**
- `--display ds9` — overlay in DS9 (DS9 must be running)
- `--display mpl` — matplotlib figure (headless-friendly)

**Analysis options:**
- `--bkg` — fit 2D polynomial background; show model + residuals
- `--phot` — run aperture photometry; overlay sources on image
- `--bkg --phot` — combined (shares same polynomial fit)

**Example: Local cache with background + photometry**

```bash
python roman_view_sca.py --display mpl --bkg --phot \
  ./cache r0003201001001001004_0001_wfi11_f106_cal.asdf
```

### `roman_phot.py`

Batch photometry across all 18 SCAs of an exposure.

**Inputs:**
- `--uri-file PATH` — text file with one S3 URI or local path per line

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

**Output:**
- `roman_phot_combined.csv` — combined source table (all SCAs)
- `roman_phot_summary.csv` — per-SCA statistics
- `roman_phot_bkg_mosaic.png` — background mosaic in focal-plane layout
- `roman_phot_histograms.png` — SNR, background level, and source count distributions
- Per-SCA CSVs (if `--per-sca`)

**Example: High-sensitivity faint source detection**

```bash
python roman_phot.py --uri-file my_exposure.txt \
  --fwhm 2.0 --det-sigma 3 --snr-threshold 2 \
  --bkg-mosaic --per-sca
```

### `roman-mast` Command

Query MAST for Roman observations.

```bash
roman-mast --program 114 --pass 57 --list
```

### `roman-fits` Command

Convert ASDF to FITS for compatibility with other tools.

```bash
roman-fits r0003201001001001004_0001_wfi11_f106_cal.asdf
```

### `roman-metadata` Command

Extract and display metadata from ASDF files.

```bash
roman-metadata r0003201001001001004_0001_wfi11_f106_cal.asdf
```

## Architecture

### Photometry Module

The core photometry functionality is in `photometry.py`, which contains the `SourcePhotometry` class:

- **Source detection** — DAOStarFinder for point source identification
- **Aperture photometry** — circular apertures with local background subtraction
- **PSF photometry** — integrated Gaussian PRF fitting
- **DS9 visualization** — region files and overlays

This replaces the external `roman-lolo` dependency, making the package self-contained.

### Streaming Pipeline

1. **Input** — S3 URI, local path, or MAST query
2. **Stream** — anonymous S3 access via s3fs, or MAST API
3. **Load** — roman_datamodels + ASDF parsing
4. **WCS** — compute SIP approximation for FITS compatibility
5. **Display** — DS9 (XPA) or matplotlib

### Focal-Plane Layout

The 18 WFI SCAs are arranged in a 3×6 grid with specific rotations:
- SCAs 3, 6, 9, 12, 15, 18 — no rotation (r=0)
- All others — 180° rotation (r=180)

Background mosaics and focal-plane visualizations respect this layout.

## Dependencies

Core:
- `numpy`, `astropy` — data structures, WCS
- `roman_datamodels`, `romancal`, `rad` — Roman-specific formats
- `photutils` — source detection, aperture photometry
- `matplotlib` — visualization
- `s3fs`, `fsspec` — S3 streaming
- `keyring`, `python-dotenv` — credential handling

Optional (for DS9):
- DS9 + XPA (local installation required)

## FAQ

**Q: Can I use this without MAST authentication?**

A: Yes! All tutorials use anonymous S3 access (`s3://stpubdata/...`). MAST authentication is only needed if querying private programs.

**Q: How do I cache files locally for fast iteration?**

A: Use `bash download_all_scas.sh` to download the tutorial dataset (~3.4 GB). Then use `--uri ./cache --filename ...` or create a local URI file with `cache/` paths.

**Q: Why are there warnings about deprecated photutils parameters?**

A: The code currently uses older photutils API (xcentroid, ycentroid, sharplo/sharphi). These will be updated in a future version to match photutils 2.0+.

**Q: Can I run this in a Docker container?**

A: Yes. Requirements: Python 3.10+, pip. No DS9 display support in containers unless running with X11 forwarding. Use `--display mpl` for headless mode.

**Q: What does the "source mask" percentage mean in background maps?**

A: Percentage of pixels masked during background fitting (bright sources, bad pixels, etc.). Values 15–20% are typical for science data.

## Contributing

Contributions welcome! For questions or issues:
1. Check existing issues on GitHub
2. Submit bug reports with example commands and output
3. Propose enhancements with motivation and use case

## License

TBD

## References

- [Roman Space Telescope Documentation](https://roman.gsfc.nasa.gov)
- [MAST Archive](https://mast.stsci.edu)
- [photutils Documentation](https://photutils.readthedocs.io)
- [astropy Documentation](https://docs.astropy.org)
