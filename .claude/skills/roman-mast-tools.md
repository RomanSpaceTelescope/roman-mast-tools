---
name: roman-mast-tools
description: Orient before working on `/home/mrizzo/roman_it_shared/roman-mast-tools`. Read this FIRST whenever a task mentions the Roman Space Telescope, WFI, MAST streaming, ASDF products, `roman_datamodels`, `astroquery.mast.MastMissions`, or edits any file in that repo. Explains which files are the clean foundation vs. the legacy layer, the two independent MAST filter axes (product_type vs kinds), the fast-path filename synthesis, the conda env + SSL cert hook, the 60-second pre-signed-URL trap, and the parallel `stream_materialized` flow. Skimming this saves you from writing code against the wrong module, triggering the slow `get_unique_product_list` path, or shipping a sink that dies on HTTP 403 halfway through.
---

# Roman MAST tools — orientation

This skill orients a fresh Claude instance on `/home/mrizzo/roman_it_shared/roman-mast-tools`. The user (Maxime Rizzo) develops tools to stream Roman WFI data from MAST for the Space Telescope's commissioning phase.

## Big picture

The user rewrote the stack **from scratch** on top of two clean modules:

- **`roman_mast.py`** — foundation. Search, filter, filename synthesis, exposure grouping, streaming. No dependency on `astropy.io.fits` or `pyds9`.
- **`roman_fits.py`** — output layer. Streams into local FITS files or into a running DS9 as a WCS mosaic. Depends on `roman_mast`.

Everything else in the repo is **legacy**:

- `query_utils.py`, `streaming_utils.py`, `peek_mast_data.py`, `write_wfi_fits.py`, `view_metadata.py`, `test_mast_adaptation.py`

**`export_metadata_csv.py` is a special case** — the CLI at the bottom is legacy, but its helpers (`flatten_metadata`, `extract_row`, `extract_rows`, `write_csv`, `write_metadata_csv`) are imported by `roman_fits.py` to drop a metadata CSV alongside every FITS/DS9 output. If you extend metadata extraction, edit the helpers at the top of that file — do NOT touch its CLI, do NOT copy the helpers elsewhere.

**Do not extend the legacy files.** If the user asks for a new tool, build it on top of `roman_mast.py` (import `list_data`, `DataResults`, `Exposure`, `stream_exposure`, `PRODUCT_KINDS`, `add_list_data_args`, `list_data_from_args`).

## Conda env

Run everything with **`conda run -n roman-mast-tools python ...`** (or activate that env first). The user also has `comm_streaming_env` — do NOT use it, it's a different environment.

**Path caveat:** conda's `conda activate` isn't on PATH by default in fresh shells here. Source it explicitly:

```bash
source /shared/spack/opt/spack/linux-zen2/miniforge3-25.3.0-3-wyxgbjkw5ewvy7ckmryhyoadzi75ufom/etc/profile.d/conda.sh
conda activate roman-mast-tools
```

**SSL cert fix — activation hook is installed.** The env's baked-in OpenSSL default cafile path (`/home/mrizzo/roman-mast-tools/ssl/cert.pem`) does not exist. `requests`-based calls (astroquery auth + `query_criteria`) work because they use `certifi` directly, but `MastMissions.read_product()` streams via `fsspec`→`aiohttp` which uses OpenSSL's default context and fails with `SSLCertVerificationError`. That failure is hidden — fsspec wraps it as a bare `FileNotFoundError(url)`, so all you see is `ERROR streaming SCA NN (...): <giant S3 pre-signed URL>` with no error class or message.

The fix is an activation hook at `/home/mrizzo/.conda/envs/roman-mast-tools/etc/conda/activate.d/ssl_cert.sh` that exports `SSL_CERT_FILE=$CONDA_PREFIX/lib/python3.12/site-packages/certifi/cacert.pem`. If streaming ever regresses to URL-only error messages, that's the first thing to check — `echo $SSL_CERT_FILE` inside an activated env; if empty, the hook is missing.

## Reference notebook

`/home/mrizzo/comm_streaming_example.ipynb` is the canonical example of how MAST streaming works. Anything you build here should reproduce its idiom cleanly:

```python
missions = MastMissions(mission='roman')
missions.login(token=...)
results  = missions.query_criteria(program=114, pass=57, detector='WFI04', select_cols=...)
products = missions.get_unique_product_list(results)   # ← SLOW; the foundation avoids this
filtered = missions.filter_products(products, file_suffix='_cal')
af       = missions.read_product(filtered['filename'][3])   # streaming reader
```

## Two independent filter axes — do not fuse them

This is the single most important thing to internalize. MAST returns two families of rows, and each row can have several file kinds. The old `data_level` argument fused them and hid mosaics; the new API separates them:

| Axis | Lever | Where filtered |
|---|---|---|
| **Which rows** MAST returns | `product_type` (`'l2'` = per-SCA, `'p_visit_coadd'` = mosaic tile) | Server-side, in `query_criteria` |
| **Which files per row** | `kinds` (list of keys from `PRODUCT_KINDS`) | Client-side, in `_synthesize_products` |

**`PRODUCT_KINDS` is the single source of truth** for the file suffixes Roman writes. Adding a new suffix is a one-line entry there and every path (fast synthesis, `--kinds` selection, `--list-kinds` help, `enumerate_products=True`) picks it up automatically. Kinds are **family-aware** — a `'cal'` kind only expands `product_type='l2'` rows, a `'coadd'` kind only expands `product_type='p_visit_coadd'` rows. Mixing them in one call (`kinds=['cal', 'coadd']`) is safe.

`data_level` is a **backwards-compat shortcut** (`1→'uncal'`, `2→'cal'`, `'gw'→'gw'`, `None→all`). If a task needs anything else — mosaics, catalogs, segmentation maps — use `kinds`, not `data_level`.

## The fast path is not optional

`missions.get_unique_product_list(results)` batches to 1000 rows per HTTP call but the server does per-row product enumeration internally. For wide searches (hundreds of rows) it will appear to hang.

**Default path** (`enumerate_products=False`): derive filenames locally from the `fileSetName` column. Roman filenames follow strict patterns:

- **Per-SCA** (`product_type='l2'`): `r{19-digit visit_id}_{4-digit exp}_wfi{SCA:02}_{filter}{suffix}{ext}`
- **Mosaic** (`product_type='p_visit_coadd'`): `r{5-digit program}_p_v{...}p{...}x{X}y{Y}_{filter}{suffix}{ext}`

Two regexes in `roman_mast.py` encode this: `_FILESET_RE` (per-SCA, in `_group_exposures`) and `_COADD_FILENAME_RE` (in `_classify_products`, disambiguates `_segm`/`_cat` between families).

Only use `enumerate_products=True` when the user is auditing for **unknown** suffixes we haven't registered — otherwise the fast path is always correct and always faster.

## visit_id decodes as PPPPPCCAAASSSOOOVVV

19 digits, packed 5/2/3/3/3/3:

- **PPPPP** program (5)
- **CC** execution_plan (2)
- **AAA** pass (3)
- **SSS** segment (3)
- **OOO** observation (3)
- **VVV** visit (3)

`parse_visit_id()` is the public helper. All six are:
1. First-class **filters** on `list_data()` and the CLI (`--program --execution-plan --pass --segment --observation --visit`). MAST accepts them as native query fields — no client-side wildcard synthesis.
2. Auto-decoded **fields** on `Exposure` (via `__post_init__`).
3. Chunk-colored **columns** in `print_summary`.

## Exposure is the display unit

`DataResults.exposures` groups rows by `(visit_id, exposure_number)` — the natural unit (18 SCAs per exposure nominally). Each `Exposure` carries:

- `visit_id`, `exposure` — identity
- `program`, `execution_plan`, `pass_`, `segment`, `observation`, `visit` — decoded components
- `optical_element`, `exposure_start_time`, `exposure_time` — display metadata
- `scas`, `filenames` — per-SCA files at the current `data_level`
- `.n_scas`, `.missing_scas` — health check (nominal 18)

Selection: `res.select(1)` (1-based index, matches `print_summary` numbering), `res.select(('visit_id', exp))`, or pass an `Exposure` through unchanged.

## DataResults carries the auth'd session

`res.missions` is the still-logged-in `MastMissions` — reuse it for streaming, don't re-auth. Convenience methods:

- `res.stream(key)` → `{sca: AsdfFile}` (uses `stream_exposure`)
- `res.to_fits(key, out_dir=..., compress=..., sip_degree=4)` — stream + write, closes buffers
- `res.to_ds9(key, dq_overlay=True)` — stream + XPA pipe, closes buffers

`close_streams(af_dict)` is the manual cleanup. All streams are `AsdfFile` objects; wrap with `rdm.open(af)` for a Roman datamodel.

These library methods are **single-consumer, single-URL-lifetime**. If you're writing library code that needs to feed the stream to multiple sinks (FITS + DS9 + CSV, etc.), use `roman_fits.stream_materialized(exp, res.missions, ..., max_workers=8)` instead — it fans out over a thread pool, materializes arrays inline, and returns a `{sca: DataModel}` dict that stays valid past URL expiry. See the streaming trap section below for why this matters.

## CLI arg block is shared

`roman_mast.add_list_data_args(parser)` + `list_data_from_args(args)` are the shared entry points. Any new CLI (`roman_export.py`, etc.) plugs in with two lines and inherits every filter for free. `roman_fits.py` already works this way.

## Diagnostics + color

- `_log(msg)` prints `[roman_mast] msg` to **stderr** (leaves stdout clean for piping).
- `_timed(label)` is a context manager that logs start/end/elapsed — wrap slow calls in it.
- `VERBOSE` module global (toggled by `--quiet`) silences diagnostics.
- `COLOR` module global (`--color always/never/auto`, respects `NO_COLOR` env var). Auto-disabled when stdout isn't a TTY. `_ansi(text, *codes)` is the low-level helper.
- Table color uses **change-highlighting**: repeats dim, changed values bold. `_CHUNK_COLORS` maps the six visit_id chunks to consistent colors used in the header row AND when coloring the raw visit_id string.

## Streaming traps to know about

Three problems bit us in real use — the fixes are in place, but the traps are structural to how MAST + fsspec + roman_datamodels compose, so any new sink can trip them again.

### 1. Pre-signed S3 URLs expire in 60 seconds

MAST-issued URLs carry `X-Amz-Expires=60`. Every `read_product(filename)` mints a fresh URL, and every byte read pulls through fsspec's HTTP cache — but the cache is keyed on byte ranges, so any block that hasn't been read yet needs a live URL to fetch. Concretely: if you open all 18 SCAs at t=0, then walk them at t=90s asking for `dm.dq[...]`, SCA 1's DQ block hasn't been read yet, and its URL is dead → HTTP 403 `Forbidden`.

**Rule:** if a stream will be consumed by more than one pass, or the passes span > ~30 s, materialize every array you'll need **into numpy** while the URL is fresh.

`roman_fits.stream_materialized(exp, missions, scas=..., max_workers=8)` does this. Per SCA: `read_product` → `rdm.open` → `_materialize_dm` → done. `_MATERIALIZE_ATTRS` at the top of `roman_fits.py` is the single-source list of what gets pulled — currently `('data', 'dq')`. Add `err` / `var_poisson` / `var_rnoise` / `var_flat` there if a future sink needs them, but each one is another 17-67 MB × 18 SCAs, so keep the list minimal.

### 2. Datamodels close their underlying AsdfFile on gc

`rdm.open(af)` returns a datamodel that holds a strong reference to `af` — but the moment the datamodel is dropped, the AsdfFile closes. If you build a dict of raw AsdfFiles and hand it to two sinks that each internally call `rdm.open(af)`, the first sink's datamodel closes the AsdfFile when it goes out of scope, and the second sink then sees "Cannot access data from closed ASDF file".

**Rule:** if a stream will be consumed by more than one sink, wrap **once** into datamodels up front and pass the `{sca: DataModel}` dict around. All sink `_open_dm` helpers already pass dm-like objects through unchanged. `stream_materialized` returns exactly this dict.

### 3. Errors from `read_product` are opaque

fsspec wraps any HTTP/SSL/etc. failure as a bare `FileNotFoundError(url)`. The exception `str()` is just the URL — no error class, no cause message. If you see `ERROR streaming SCA NN (...): https://s3...`, that is the error being surfaced, not a truncated log line. Real causes we've seen:

- **SSL cert path broken** → see the conda-env section above.
- **URL expired** → the 60-second rule above.
- **AWS token revoked mid-stream** → same shape; retry the whole stream.

`_stream_one_sca` in `roman_fits.py` logs `{type(e).__name__}: {e}` to at least surface the exception class. When debugging, wrap `read_product` in `try/except` and print `traceback.format_exc()` — the real cause is chained under `__cause__`.

## Parallel streaming

`stream_materialized` fans out over a `ThreadPoolExecutor` (default `max_workers=8`). Threads are the right shape here because each `read_product` is a blocking network read that releases the GIL, and MAST mints an independent pre-signed URL per request, so there's nothing shared to contend on. Empirically 8 workers takes an 18-SCA exposure from ~2m40s sequential to ~16s.

- CLI flag: `--workers N` on `roman_fits.py`. Default 8, set 1 for the sequential path (still preserved verbatim in the function body — one less thing to debug when a machine misbehaves).
- Results are re-sorted into ascending SCA order before returning, so downstream logs / CSV rows come out `01…18` regardless of who finished first.
- If a machine is bandwidth-shared or hits S3 throttling, drop to `--workers 4`. Diminishing returns past ~8 in current tests.

## The `roman_fits.py` CLI is the reference multi-sink flow

`python roman_fits.py --program ... --to {fits,ds9}` is the canonical example of streaming one exposure and dispatching to multiple sinks. The current flow per exposure is:

1. `stream_materialized(exp, res.missions, scas=..., max_workers=args.workers)` — parallel fetch, arrays materialized inline. Returns `{sca: DataModel}`.
2. `write_metadata_csv(dm_dict, exp, output=...)` — unless `--no-metadata`. Runs first because it's cheap and self-contained.
3. `to_fits_files(dm_dict, exp, ...)` OR `to_ds9(dm_dict, exp, ...)` — image sink.
4. `close_streams(dm_dict)` in a `finally`.

If you add a third sink, follow the same shape: take a `dm_dict` (not an af_dict, not a res+key pair), never re-stream, don't call `rdm.open` yourself. Put opt-out flags on the CLI (`--no-metadata` is the template) so users can skip sinks they don't need.

**Why the metadata CSV runs on every FITS/DS9 output by default:** we're already streaming the whole exposure into memory — extracting metadata is essentially free (~200 ms per exposure) and the CSV is often the artifact the user actually goes back to. Filename convention: `metadata_{visit_id}_exp{NN}.csv`, dropped into `--out-dir` in fits mode or cwd in ds9 mode (overridable with `--metadata-dir`).

## What's deliberately deferred

- **Guide-star region overlay** — old `stream_to_ds9` emitted an fk5 polygon per SCA for the guide-window window. Port as `guide_star_regions(dm_dict) → region text` when needed.
- **Interactive prompts** — the old CLIs had `prompt_query_params`. The new stack replaces this with the `--list` → adjust → re-run pattern (`roman_fits.py --list`). Don't add interactive prompts back.
- **Per-SCA streaming pipeline (option: stream one SCA → drain through all sinks → discard → next).** Would cut peak RAM from ~1.5 GB (all 18 SCAs' data+dq) to ~85 MB (one SCA). Not implemented because DS9's mosaic API (`fits mosaicimage wcs` via XPA) wants the whole MEF at once — going per-SCA would force two divergent code paths (fits vs ds9). Leave until we hit a real memory ceiling.

## Working style

- The user prefers **structured filters over string parsing**. If you find yourself writing `'1-3,5'` parsing, the answer is `parse_int_spec()` in `roman_mast.py` — already exists.
- **New CLIs get a `--list` mode** that prints the summary without side effects. Encourages recon before commit.
- Progress bars: **tqdm on streaming loops** (over SCAs). Not on metadata queries.
- **Sanity-check every non-trivial addition end-to-end.** The user's data is real (program 114, pass 57 is the reference MRT-7b test set) — running a query against it takes ~5 seconds and catches many mistakes. Suggested commands after a change:
  ```bash
  conda run -n roman-mast-tools python roman_mast.py --program 114 --pass 57 --sca-only --max-rows 6
  conda run -n roman-mast-tools python roman_fits.py --program 114 --pass 57 --sca-only --list
  ```
- Style: **short comments explaining WHY**, matching the density in `roman_mast.py`. Avoid boilerplate docstrings on trivial helpers. Do not add heavyweight docstring headers to short functions.

## Repository layout at a glance

```
roman-mast-tools/
├── roman_mast.py           ← foundation. Auth, filters, list_data, exposures, sequential streaming.
├── roman_fits.py           ← output layer. stream_materialized, to_fits_files, to_ds9, CLI.
├── export_metadata_csv.py  ← MIXED: helpers (flatten_metadata, extract_row(s),
│                              write_csv, write_metadata_csv) are current and imported
│                              by roman_fits.py. Its bottom-of-file CLI is legacy.
├── .claude/skills/roman-mast-tools.md   ← this file
│
├── comm_streaming_example.ipynb   ← reference notebook (in /home/mrizzo/, not here)
│
└── (legacy — do not extend)
    query_utils.py streaming_utils.py peek_mast_data.py
    write_wfi_fits.py view_metadata.py test_mast_adaptation.py
```

## When something surprises you

- If a query returns a row family you don't recognize → check `product_type` on the metadata result. Two known values: `'l2'`, `'p_visit_coadd'`. Anything else is new and worth logging.
- If `_classify_products` puts something in `'other'` → that filename ends with a suffix not in `PRODUCT_KINDS`. Add the entry; that's the intended signal.
- If the fast path fails but `enumerate_products=True` works → `fileSetName` was missing from `results.colnames`. Check the `columns=` list.
- If `Exposure.__post_init__` leaves the decoded fields as `None` → the visit_id isn't 19 digits. Rare in real MAST data but not impossible for engineering datasets.
- **If streaming logs `ERROR streaming SCA NN (...): <full S3 URL>` with no error class** → that's fsspec masking the real cause. Check `echo $SSL_CERT_FILE` first (cert-path breakage on this machine has bitten us; activation hook should be installed). If SSL is fine, it's URL expiry — you're consuming a stream more than 60 s after opening it. Switch to `stream_materialized`.
- **If a second consumer sees `Cannot access data from closed ASDF file`** → an earlier consumer's `rdm.open(af)` datamodel went out of scope and closed the underlying AsdfFile. Wrap into `{sca: DataModel}` once, up front, and pass the dm_dict — see the streaming-traps section.
- **If streaming is slow (>10 s per SCA)** → check `--workers`. Default is 8; if someone downgraded to 1 or the machine is bandwidth-shared, single-stream is ~9 s/SCA and 8-stream is ~1 s/SCA amortized.
