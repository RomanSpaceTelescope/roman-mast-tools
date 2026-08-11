# Plan: anonymous S3 streaming (`roman_s3.py`)

## Goal

Stream Roman WFI data directly from a public S3 bucket (no MAST auth) for
offline practice and tutorial datasets. Example bucket:

```
s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/
```

That prefix holds one full 18-SCA visit (`r0003201001001001004`, 4 exposures,
filter F106, all kinds: cal/uncal/wcs/segm/cat) plus a mosaic and a handful of
single-SCA extras. Files follow the standard Roman naming convention so the
existing regex machinery applies directly.

## Architecture

New **`roman_s3.py`** module — no edits to `roman_mast.py`, `roman_fits.py`,
or any legacy file.

```
list_s3_data(s3_uri, ...)         → [Exposure]          (discovery)
stream_s3_materialized(exp, uri)  → {sca: DataModel}     (same shape as
                                                           roman_fits.stream_materialized)
     ↓
to_fits_files / to_ds9            unchanged  (take {sca: DataModel})
```

## `roman_s3.py` contents

### `list_s3_data(s3_uri, *, kinds=None, data_level=2, scas=None) → list[Exposure]`

- `s3fs.S3FileSystem(anon=True).ls(s3_uri)` to list filenames
- Parse basenames with `_FILESET_RE` (imported from `roman_mast`)
- Group into `Exposure` objects via the same logic as `_group_exposures`
- Filter by `kinds`/`data_level` using `_resolve_kinds` (imported)
- Returns a plain `list[Exposure]` — no `DataResults` wrapper (no session to carry)
- Helper `print_s3_summary(exposures)` for the `--list` recon step

### `stream_s3_materialized(exposure, s3_base_uri, *, scas=None, max_workers=8, show_progress=True) → dict[int, DataModel]`

- ThreadPoolExecutor over SCAs (same structure as `roman_fits.stream_materialized`)
- Per SCA: `s3fs.S3FileSystem(anon=True).open(s3_base_uri + filename)`
  → `asdf.open(fobj)` → `rdm.open(af)` → `_materialize_dm(dm)` (imported from `roman_fits`)
- Returns `{sca: DataModel}` — drop-in replacement for `stream_materialized` output
- No 60-second URL trap (public bucket, no pre-signed URLs), but still materialize
  eagerly — lazy S3 reads are slow

### CLI

```
python roman_s3.py --s3-uri s3://... [--list] [--to fits|ds9]
                   [--exposures N] [--scas N] [--workers N]
                   [--out-dir DIR] [--compress] ...
```

- `--s3-uri`   required; the S3 prefix ending in `/`
- `--list`     print summary and exit (recon-first pattern)
- All other flags mirror `roman_fits.py` CLI; calls the same
  `to_fits_files` / `to_ds9` / `write_metadata_csv` sinks

## What is reused unchanged

| Symbol | From |
|---|---|
| `Exposure`, `_FILESET_RE`, `_group_exposures` | `roman_mast` |
| `PRODUCT_KINDS`, `_resolve_kinds` | `roman_mast` |
| `parse_int_spec`, `add_list_data_args` (if shared CLI flags desired) | `roman_mast` |
| `_materialize_dm` | `roman_fits` |
| `to_fits_files`, `to_ds9` | `roman_fits` |
| `write_metadata_csv` | `export_metadata_csv` |

## Sanity-check command

```bash
conda run -n roman-mast-tools python roman_s3.py \
  --s3-uri s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
  --list
```

Then to stream one exposure to FITS:

```bash
conda run -n roman-mast-tools python roman_s3.py \
  --s3-uri s3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.2/ \
  --exposures 1 --to fits --out-dir /tmp/s3_test
```

## Key differences vs the MAST path

| | MAST | S3 |
|---|---|---|
| Auth | MAST token required | anonymous |
| Discovery | `query_criteria` (server-side) | `fs.ls(uri)` (client-side) |
| File open | `missions.read_product(filename)` | `fs.open(s3://...)` |
| URL expiry | 60-second pre-signed URLs → must materialize | n/a, but still materialize for speed |
| Session object | `MastMissions` carried in `DataResults` | `s3_base_uri` string passed to stream fn |
