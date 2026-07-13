---
name: roman-mast-tools
description: Orient before working on `/home/mrizzo/roman_it_shared/roman-mast-tools`. Read this FIRST whenever a task mentions the Roman Space Telescope, WFI, MAST streaming, ASDF products, `roman_datamodels`, `astroquery.mast.MastMissions`, or edits any file in that repo. Explains which files are the clean foundation vs. the legacy layer, the two independent MAST filter axes (product_type vs kinds), the fast-path filename synthesis, and the conda env. Skimming this saves you from writing code against the wrong module or triggering the slow `get_unique_product_list` path.
---

# Roman MAST tools — orientation

This skill orients a fresh Claude instance on `/home/mrizzo/roman_it_shared/roman-mast-tools`. The user (Maxime Rizzo) develops tools to stream Roman WFI data from MAST for the Space Telescope's commissioning phase.

## Big picture

The user rewrote the stack **from scratch** on top of two clean modules:

- **`roman_mast.py`** — foundation. Search, filter, filename synthesis, exposure grouping, streaming. No dependency on `astropy.io.fits` or `pyds9`.
- **`roman_fits.py`** — output layer. Streams into local FITS files or into a running DS9 as a WCS mosaic. Depends on `roman_mast`.

Everything else in the repo is **legacy**:

- `query_utils.py`, `streaming_utils.py`, `peek_mast_data.py`, `write_wfi_fits.py`, `export_metadata_csv.py`, `view_metadata.py`, `test_mast_adaptation.py`

**Do not extend the legacy files.** If the user asks for a new tool, build it on top of `roman_mast.py` (import `list_data`, `DataResults`, `Exposure`, `stream_exposure`, `PRODUCT_KINDS`, `add_list_data_args`, `list_data_from_args`).

## Conda env

Run everything with **`conda run -n roman-mast-tools python ...`** (or activate that env first). The user also has `comm_streaming_env` — do NOT use it, it's a different environment.

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

## CLI arg block is shared

`roman_mast.add_list_data_args(parser)` + `list_data_from_args(args)` are the shared entry points. Any new CLI (`roman_export.py`, etc.) plugs in with two lines and inherits every filter for free. `roman_fits.py` already works this way.

## Diagnostics + color

- `_log(msg)` prints `[roman_mast] msg` to **stderr** (leaves stdout clean for piping).
- `_timed(label)` is a context manager that logs start/end/elapsed — wrap slow calls in it.
- `VERBOSE` module global (toggled by `--quiet`) silences diagnostics.
- `COLOR` module global (`--color always/never/auto`, respects `NO_COLOR` env var). Auto-disabled when stdout isn't a TTY. `_ansi(text, *codes)` is the low-level helper.
- Table color uses **change-highlighting**: repeats dim, changed values bold. `_CHUNK_COLORS` maps the six visit_id chunks to consistent colors used in the header row AND when coloring the raw visit_id string.

## What's deliberately deferred

- **CSV metadata sidecar** — the old `write_wfi_fits.py` wrote a per-exposure metadata CSV via `flatten_metadata` from `export_metadata_csv.py`. The user said "ask later". If they do, add it as `roman_fits.write_metadata_csv(af_dict, out_dir, exposure)` — same style as the other output helpers. The old `flatten_metadata` and its priority-key ordering are worth copying.
- **Guide-star region overlay** — old `stream_to_ds9` emitted an fk5 polygon per SCA for the guide-window window. Port as `guide_star_regions(af_dict) → region text` when needed.
- **Interactive prompts** — the old CLIs had `prompt_query_params`. The new stack replaces this with the `--list` → adjust → re-run pattern (`roman_fits.py --list`). Don't add interactive prompts back.

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
├── roman_mast.py       ← foundation. Auth, filters, list_data, exposures, streaming.
├── roman_fits.py       ← output layer. to_fits_files, to_ds9. CLI.
├── .claude/skills/roman-mast-tools.md   ← this file
│
├── comm_streaming_example.ipynb   ← reference notebook (in /home/mrizzo/, not here)
│
└── (legacy — do not extend)
    query_utils.py streaming_utils.py peek_mast_data.py
    write_wfi_fits.py export_metadata_csv.py
    view_metadata.py test_mast_adaptation.py
```

## When something surprises you

- If a query returns a row family you don't recognize → check `product_type` on the metadata result. Two known values: `'l2'`, `'p_visit_coadd'`. Anything else is new and worth logging.
- If `_classify_products` puts something in `'other'` → that filename ends with a suffix not in `PRODUCT_KINDS`. Add the entry; that's the intended signal.
- If the fast path fails but `enumerate_products=True` works → `fileSetName` was missing from `results.colnames`. Check the `columns=` list.
- If `Exposure.__post_init__` leaves the decoded fields as `None` → the visit_id isn't 19 digits. Rare in real MAST data but not impossible for engineering datasets.
