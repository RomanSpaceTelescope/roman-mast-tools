# Installation

These instructions create a Python environment suitable for running the
scripts in this repo and the `comm_streaming_example.ipynb` notebook.

The environment installs `roman_datamodels`, `rad`, `fsspec[s3]`,
`matplotlib`, and a pre-release build of `astroquery` from a specific pull
request that adds the Roman MAST search/streaming support the tools depend
on. `romancal` is intentionally not installed — nothing in this repo imports
it, and it depends on `romanisim` → `galsim`, which cannot be built on
Windows.

## Prerequisites

- Python 3.12 (any recent 3.10+ should also work; 3.12 is what we test with)
- `git` on your `PATH` (pip needs it to clone the astroquery PR)
- A C toolchain in case any dependency has to build from source (`gcc`,
  `make`; already present on the shared cluster)

## Quick start (conda / mamba) — recommended

The simplest approach. An `environment.yml` is included with all dependencies:

```bash
mamba env create -f environment.yml
mamba activate roman-mast-tools

# Optional, for the notebook:
pip install jupyterlab ipykernel
python -m ipykernel install --user \
  --name roman-mast-tools \
  --display-name "Python (roman-mast-tools)"
```

Then verify the install:

```bash
python -c "
import roman_datamodels, astroquery, fsspec, matplotlib
from astroquery.mast import MastMissions
print('roman_datamodels', roman_datamodels.__version__)
print('astroquery      ', astroquery.__version__)
print('OK')
"
```

## Windows (miniforge)

`environment.yml` works on Windows as-is once `romancal` is out of the
picture (see above). Note that `astroquery` used to arrive transitively via
`romancal`; it is now an explicit dependency in `pyproject.toml`, so
`pip install -e .` pulls the pinned pre-release build on every platform.

If you are on an older checkout that still lists `romancal`, or
`pip install -e .` tries to build `galsim` and fails with
`ValueError: list.remove(x): x not in list`, install the deps explicitly and
skip dependency resolution for the editable install:

```powershell
mamba create -n roman-mast-tools python=3.12 numpy matplotlib astropy tqdm keyring python-dotenv requests s3fs git pip
mamba activate roman-mast-tools
pip install "roman_datamodels>=1.0.0" "rad>=1.0.0" "fsspec[s3]" photutils reproject scikit-image astroalign "git+https://github.com/snbianco/astroquery.git@2026.2"
pip install --no-deps -e .
```

DS9/XPA is awkward on Windows; use `--to fits` and `--display mpl` there.

## Alternative: Step-by-step (venv + pip)

This is the recommended path. It uses the standard-library `venv` module
and pip; no conda required.

1. Clone the repository (if you haven't already) and change into it:

   ```bash
   git clone <repo-url> roman-mast-tools
   cd roman-mast-tools
   ```

2. Create the virtual environment. Pick any location you like; here we put
   it next to the repo so it is easy to find but is not committed:

   ```bash
   python3 -m venv ~/roman-mast-tools-env
   ```

   On Windows use `py -3.12 -m venv %USERPROFILE%\roman-mast-tools-env`.

3. Activate the environment:

   ```bash
   source ~/roman-mast-tools-env/bin/activate
   ```

   On Windows: `%USERPROFILE%\roman-mast-tools-env\Scripts\activate` (cmd)
   or `& $env:USERPROFILE\roman-mast-tools-env\Scripts\Activate.ps1`
   (PowerShell).

4. Upgrade the packaging tools inside the environment (avoids surprises
   when building `astroquery` from the PR branch):

   ```bash
   pip install --upgrade pip setuptools wheel
   ```

5. Install the package and its dependencies (editable, so edits to the
   scripts take effect without reinstalling):

   ```bash
   pip install -e .
   ```

   The `astroquery @ git+https://...` dependency in `pyproject.toml` will
   build `astroquery` from the pinned pre-release branch. This can take a
   minute or two.

6. (Optional, only if you want to run the notebook) Install JupyterLab and
   register a kernel for this environment:

   ```bash
   pip install jupyterlab ipykernel
   python -m ipykernel install --user \
     --name roman-mast-tools \
     --display-name "Python (roman-mast-tools)"
   ```

7. Verify the install:

   ```bash
   python -c "
   import roman_datamodels, astroquery, fsspec, matplotlib
   from astroquery.mast import MastMissions
   print('roman_datamodels', roman_datamodels.__version__)
   print('astroquery      ', astroquery.__version__)
   print('OK')
   "
   ```

   The `astroquery` version will be a dev string like
   `0.4.12.dev583+ga5214b2ba` — that confirms the PR branch was installed
   rather than a release from PyPI.

## Running the notebook

With the environment activated:

```bash
jupyter lab
```

Open `comm_streaming_example.ipynb` and choose the
**Python (roman-mast-tools)** kernel from the kernel picker.

The notebook needs a MAST auth token. Set one of:

- Environment variable (recommended, set *before* launching Jupyter):

  ```bash
  export MAST_API_TOKEN=<your-token>
  ```

- Or a file named `mast_api_token.txt` in the same directory as the
  notebook containing just the token.

## Running the scripts

The scripts in this repo (e.g. `write_wfi_fits.py`, `peek_mast_data.py`,
`display_mosaic_wfi.py`) can be run directly once the environment is
active:

```bash
source ~/roman-mast-tools-env/bin/activate
python write_wfi_fits.py --help
```

## Updating

The astroquery PR is expected to be merged and released before launch and
commissioning. Once that happens, the `astroquery @ git+…` line in
`pyproject.toml` (and `environment.yml`) can be replaced with a normal
`astroquery>=<version>` pin.

To reinstall from scratch:

```bash
deactivate                            # if the env is active
rm -rf ~/roman-mast-tools-env
python3 -m venv ~/roman-mast-tools-env
source ~/roman-mast-tools-env/bin/activate
pip install --upgrade pip setuptools wheel
pip install -e .
```

## Troubleshooting

- **`error: could not find git`** — install git and re-run step 5.
- **`ModuleNotFoundError: astroquery.mast`** after install — either you're
  on the system Python, not the venv (re-run the activate command from step
  3 and check `which python` / `where python`), or you installed with
  `--no-deps` on an old checkout and need to `pip install
  "git+https://github.com/snbianco/astroquery.git@2026.2"` by hand.
- **Notebook kernel missing** — re-run step 6, then reload the JupyterLab
  browser tab.
- **`MAST token not found!`** — see the notebook section above.
