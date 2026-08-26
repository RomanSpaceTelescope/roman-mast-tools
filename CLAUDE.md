Only read files I specify, and if you feel you need to read another file ALWAYS ask first

## Orientation skill and command

- `.claude/skills/roman-mast-tools.md` — detailed orientation skill covering `roman_mast.py`/`roman_fits.py`/`roman_metadata.py`/`roman_view_sca.py`/`roman_phot.py` architecture, the two MAST filter axes, streaming traps (60s pre-signed URL expiry, AsdfFile gc-on-close), conda env + SSL cert hook, and working conventions. Load it before working on anything in that module set — it does NOT currently cover `roman_telem.py`/`roman_telem_plot.py`.
- `.claude/commands/roman-mast-tools.md` — slash command that loads the above skill for orientation.

## File map

- `roman_telem.py` — CLI/API for querying Roman telemetry mnemonics from the MAST Engineering DB. Contains an inlined `MASTEngDBQuery` class (query + retry + save logic) plus YAML "plot groups" handling and the `roman-telem` CLI entry point.
- `roman_telem_plot.py` — plotting helpers (`plot_telemetry`, `plot_mnemonics`, `plot_grouped`) that consume long-format DataFrames (one row per mnemonic/time/value) with a `mnemonic` column.
- `mast_eng_db_query.py` — original standalone module `roman_telem.py` was based on; largely superseded by the inlined copy but still used elsewhere (check before assuming it's dead).
- `roman_mast.py`, `roman_fits.py`, `roman_metadata.py`, `roman_view_sca.py`, `roman_phot.py`, `photometry.py` — MAST/FITS access, metadata, SCA visualization, and photometry tools (separate CLI entry points, see `pyproject.toml`).
- `*.yaml` (e.g. `tlm_groups.yaml`) — telemetry "plot groups" config files (see schema below).

## Telemetry plot-groups YAML schema

```yaml
groups:
  <group_name>:
    label: "<display label>"
    y_label: "<axis label, optional>"
    mnemonics:
      - MNEMONIC_1
      - MNEMONIC_2
```

Used by `--plot-groups`/`--select-groups` in `roman-telem` to filter and label subplots.

## Dev workflow

- After adding/removing a top-level `.py` module used as a script or import target, update `py-modules` in `pyproject.toml`, then run `pip install -e .` to pick up the change (editable installs don't auto-discover new modules).

## MAST EDB constraints

- `--server mast` (default) hits `mast.stsci.edu`; `--server int` hits `mastint.stsci.edu` (internal/test server).
- API tokens come from `MAST_API_TOKEN` (mast) / `MAST_API_TOKEN_INT` (int) env vars, or `--api-token`.