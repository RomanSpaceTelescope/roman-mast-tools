"""Roman MAST tools — clean foundation for streaming Roman WFI data.

This module is the starting point for the rewritten roman-mast-tools stack.
Every future tool (streaming, plotting, export, ...) should build on the
primitives here rather than re-implementing them.

Public surface
--------------
    connect(token=None)          → authenticated MastMissions session
    list_data(**filters)         → search MAST, return DataResults
    print_summary(res)           → human-readable dump of a DataResults
    stream_exposure(exp, ...)    → stream every SCA of one exposure into memory
    DataResults                  → dataclass with results / products / filenames
                                    / exposures  (grouped per (visit_id, exp))
                                    + .select(key) / .stream(key) helpers
                                    + .to_fits(key, ...) / .to_ds9(key, ...)
                                      one-liners (see roman_fits.py)
    Exposure                     → one exposure: visit_id, exp #, filter, SCAs,
                                    filenames — the natural display unit

Every filter argument is optional. Omitting a filter means "return whatever
MAST has for that field" — passing no filters at all returns every Roman
product MAST knows about (be prepared for it to be large / slow).

Filter reference
----------------
    program         : int             — APT program ID (e.g. 114)
    execution_plan  : int             — execution plan within the program
    pass_           : int             — pass within the execution plan (e.g. 57)
    segment         : int             — segment within the pass
    observation     : int             — observation within the segment
    visit           : int             — visit within the observation
    detector        : str | int       — 'WFI04' / 'wfi04' / 4  → 'WFI04'
    visit_id        : str             — full ID or wildcard, e.g. '0011401057*'
    exposure        : int             — matches last 4 digits of observation_id
    optical_element : str             — e.g. 'F062', 'F129'
    exposure_type   : str             — e.g. 'WFI_IMAGE', 'WFI_DARK'
    product_type    : str             — 'l2' (raw per-SCA exposures) or
                                        'p_visit_coadd' (mosaic tiles). Set
                                        via sca_only=True as a shortcut.
    sca_only        : bool            — shortcut for product_type='l2'; drops
                                        mosaic tiles from the results.
    data_level      : 1 | 2 | 'gw' | None
                                        1 → _uncal, 2 → _cal, 'gw' → _gw,
                                        None → no product-suffix filter (returns
                                        every product kind: cal, uncal, cat,
                                        wcs, segm, ...)

Example — reproducing the comm_streaming_example.ipynb search
-------------------------------------------------------------
    >>> from roman_mast import list_data, print_summary
    >>> res = list_data(program=114, pass_=57, detector='WFI04')
    >>> print_summary(res)
    >>> res.filenames[:3]
    ['r0011401057001001001_0001_wfi04_f062_cal.asdf', ...]
"""

# Null keyring so headless envs don't hit DBus/SecretService.
import keyring
import keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

import os
import re
import sys
import time
import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Optional

from astroquery.mast import MastMissions


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

# Toggle at runtime with roman_mast.VERBOSE = False (or --quiet on the CLI).
VERBOSE = True


def _log(msg):
    if VERBOSE:
        print(f"[roman_mast] {msg}", file=sys.stderr, flush=True)


@contextmanager
def _timed(label):
    """Context manager: log start/end + elapsed for a MAST call."""
    _log(f"{label} ...")
    t0 = time.monotonic()
    try:
        yield
    finally:
        _log(f"{label} done in {time.monotonic() - t0:.1f}s")


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Columns requested from MAST. Matches the notebook plus a couple extras
# (visit_id, observation_id) that downstream tools need for grouping.
DEFAULT_COLUMNS = [
    'fileSetName',
    'program', 'execution_plan', 'pass', 'segment',
    'visit', 'observation', 'observation_id', 'visit_id',
    'optical_element', 'exposure_type', 'instrument_name', 'detector',
    'productLevel', 'product_type',
    'exposure_time', 'exposure_start_time', 'exposure_end_time',
]

# Every Roman WFI product kind we know how to build a filename for. A "kind"
# is one row in this registry; the key is the short name callers use.
#
#   suffix        — appended to fileSetName (e.g. '_cal' → 'r..._cal.asdf')
#   ext           — file extension including the dot
#   product_type  — which MAST `product_type` row family produces this file.
#                   `l2` = per-SCA exposures; `p_visit_coadd` = mosaic tiles.
#                   Kinds only get synthesized for rows of the matching family,
#                   so we never invent a `_cal.asdf` for a coadd row.
#   description   — one-line human summary for --list-kinds and help text.
PRODUCT_KINDS = {
    # --- per-SCA products (product_type='l2') ---
    'uncal': {'suffix': '_uncal', 'ext': '.asdf',    'product_type': 'l2',
              'description': 'L1 uncalibrated ramp'},
    'cal':   {'suffix': '_cal',   'ext': '.asdf',    'product_type': 'l2',
              'description': 'L2 calibrated rate image'},
    'wcs':   {'suffix': '_wcs',   'ext': '.asdf',    'product_type': 'l2',
              'description': 'L2 WCS solution'},
    'gw':    {'suffix': '_gw',    'ext': '.asdf',    'product_type': 'l2',
              'description': 'Guide-window data'},
    # --- L4 per-SCA overlays (still produced by product_type='l2') ---
    'segm_sca': {'suffix': '_segm', 'ext': '.asdf',    'product_type': 'l2',
                 'description': 'L4 segmentation map (per-SCA)'},
    'cat_sca':  {'suffix': '_cat',  'ext': '.parquet', 'product_type': 'l2',
                 'description': 'L4 source catalog (per-SCA)'},
    # --- Mosaic / coadd products (product_type='p_visit_coadd') ---
    'coadd': {'suffix': '_coadd', 'ext': '.asdf',    'product_type': 'p_visit_coadd',
              'description': 'L3 visit-coadd mosaic tile'},
    'segm':  {'suffix': '_segm',  'ext': '.asdf',    'product_type': 'p_visit_coadd',
              'description': 'L4 segmentation map (coadd)'},
    'cat':   {'suffix': '_cat',   'ext': '.parquet', 'product_type': 'p_visit_coadd',
              'description': 'L4 source catalog (coadd)'},
}

# Backward-compat shortcut: the old data_level=1/2/'gw' → one PRODUCT_KINDS key.
DATA_LEVEL_KIND = {1: 'uncal', 2: 'cal', 'gw': 'gw'}

# Legacy exports — a few callers still reference these; keep for compatibility.
DATA_LEVEL_SUFFIX = {lvl: PRODUCT_KINDS[k]['suffix']
                     for lvl, k in DATA_LEVEL_KIND.items()}
DATA_LEVEL_FILE = {lvl: (PRODUCT_KINDS[k]['suffix'], PRODUCT_KINDS[k]['ext'])
                   for lvl, k in DATA_LEVEL_KIND.items()}


def _resolve_kinds(kinds, data_level):
    """Normalize the (kinds, data_level) arguments into a list of kind names.

    Precedence:
      * kinds='all'  → every kind in PRODUCT_KINDS
      * kinds set    → validated list of names from PRODUCT_KINDS
      * data_level ∈ {1, 2, 'gw'}  → the single legacy kind
      * data_level=None + kinds=None  → all kinds (i.e. every product)
    """
    if kinds is not None:
        if isinstance(kinds, str):
            if kinds.lower() in ('all', '*'):
                return list(PRODUCT_KINDS)
            kinds = [kinds]
        kinds = list(kinds)
        bad = [k for k in kinds if k not in PRODUCT_KINDS]
        if bad:
            raise ValueError(
                f"Unknown product kind(s) {bad}. "
                f"Valid kinds: {sorted(PRODUCT_KINDS)}"
            )
        return kinds

    if data_level is None:
        return list(PRODUCT_KINDS)

    if data_level in DATA_LEVEL_KIND:
        return [DATA_LEVEL_KIND[data_level]]

    raise ValueError(
        f"data_level must be one of {list(DATA_LEVEL_KIND)} or None "
        f"(got {data_level!r}). Use the `kinds` argument for other products "
        f"({sorted(PRODUCT_KINDS)})."
    )


# ---------------------------------------------------------------------------
# Auth
# ---------------------------------------------------------------------------

def get_token(token_file='mast_api_token.txt'):
    """Return a MAST token from $MAST_API_TOKEN or *token_file* (or None)."""
    token = os.getenv('MAST_API_TOKEN')
    if token:
        return token.strip() or None
    try:
        with open(token_file) as f:
            return f.read().strip() or None
    except FileNotFoundError:
        return None


def connect(token=None):
    """Return an authenticated MastMissions session for mission='roman'."""
    if token is None:
        token = get_token()
    if not token:
        raise RuntimeError(
            "MAST token not found. Set MAST_API_TOKEN or place your token in "
            "mast_api_token.txt (see https://auth.mast.stsci.edu/info)."
        )
    missions = MastMissions(mission='roman')
    missions.login(token=token)
    return missions


# ---------------------------------------------------------------------------
# Filter normalization
# ---------------------------------------------------------------------------

def _norm_detector(det):
    """Accept 'WFI04' / 'wfi04' / 4 / '4' → 'WFI04'."""
    if det is None or det == '':
        return None
    s = str(det).strip().upper()
    if s.startswith('WFI'):
        return s
    return f'WFI{int(s):02d}'


def _norm_exposure(exp):
    """Accept int or str → MAST wildcard '*NNNN' on observation_id."""
    if exp is None or exp == '':
        return None
    return f'*{int(exp):04d}'


# kwarg → (MAST field name, normalizer). None normalizer means passthrough.
# The visit_id-component fields (execution_plan / pass / segment / observation /
# visit) all narrow server-side — MAST decodes the visit_id itself. `pass_` and
# `execution_plan` have trailing/embedded underscores because Python.
_FILTER_MAP = {
    'program':         ('program',         None),
    'execution_plan':  ('execution_plan',  None),
    'pass_':           ('pass',            None),
    'segment':         ('segment',         None),
    'observation':     ('observation',     None),
    'visit':           ('visit',           None),
    'detector':        ('detector',        _norm_detector),
    'visit_id':        ('visit_id',        None),
    'exposure':        ('observation_id',  _norm_exposure),
    'optical_element': ('optical_element', None),
    'exposure_type':   ('exposure_type',   None),
    'product_type':    ('product_type',    None),
}


def _build_search(**kwargs):
    """Turn kwargs into a MAST search dict, dropping None/blank values."""
    out = {}
    for key, val in kwargs.items():
        if val is None or val == '':
            continue
        if key not in _FILTER_MAP:
            raise TypeError(f"Unknown filter: {key!r}")
        mast_key, norm = _FILTER_MAP[key]
        out[mast_key] = norm(val) if norm else val
    return out


# ---------------------------------------------------------------------------
# Results container
# ---------------------------------------------------------------------------

# fileSetName is r{19-digit visit_id}_{4-digit exp}_wfi{SCA}_{filter}
_FILESET_RE = re.compile(
    r'r(?P<visit_id>\d{19})_(?P<exposure>\d{4})_wfi(?P<sca>\d{2})_(?P<filter>\w+)'
)

# Roman visit_id = PPPPPCCAAASSSOOOVVV (19 digits)
#   PPPPP program (5)   CC execution plan (2)   AAA pass (3)
#   SSS segment (3)     OOO observation (3)     VVV visit (3)
_VISIT_ID_RE = re.compile(
    r'^(?P<program>\d{5})(?P<execution_plan>\d{2})(?P<pass_>\d{3})'
    r'(?P<segment>\d{3})(?P<observation>\d{3})(?P<visit>\d{3})$'
)


def parse_visit_id(visit_id):
    """Split a 19-digit Roman visit_id into its six numeric components.

    Returns a dict: program, execution_plan, pass, segment, observation, visit
    (all ints). Raises ValueError on a malformed input.
    """
    m = _VISIT_ID_RE.match(str(visit_id))
    if not m:
        raise ValueError(
            f"visit_id {visit_id!r} is not 19 digits (PPPPPCCAAASSSOOOVVV)"
        )
    d = m.groupdict()
    # 'pass' is a Python keyword — normalize to 'pass' in the returned dict
    # (callers can `d['pass']`), but the regex uses 'pass_'.
    d['pass'] = d.pop('pass_')
    return {k: int(v) for k, v in d.items()}


@dataclass
class Exposure:
    """One Roman exposure — the unit you want to display / iterate over.

    An exposure has up to 18 SCAs (WFI01–WFI18). All SCAs of one exposure
    share visit_id + exposure_number + filter + start time.

    The `visit_id` is decoded into its six numeric fields (program,
    execution_plan, pass, segment, observation, visit) per the Roman
    convention PPPPPCCAAASSSOOOVVV — populated in __post_init__.
    """
    visit_id: str
    exposure: int                       # 1-based exposure number in the visit
    optical_element: Optional[str]      # e.g. 'F062'
    exposure_start_time: Any = None
    exposure_time: Any = None
    scas: list = field(default_factory=list)         # sorted list of SCA ints
    filenames: list = field(default_factory=list)    # filenames for this exposure

    # visit_id components — decoded once from the 19-digit visit_id.
    program: Optional[int]        = None
    execution_plan: Optional[int] = None
    pass_: Optional[int]          = None   # trailing underscore: `pass` is reserved
    segment: Optional[int]        = None
    observation: Optional[int]    = None
    visit: Optional[int]          = None

    def __post_init__(self):
        try:
            parts = parse_visit_id(self.visit_id)
        except ValueError:
            return   # non-standard visit_id → leave the fields as None
        self.program        = parts['program']
        self.execution_plan = parts['execution_plan']
        self.pass_          = parts['pass']
        self.segment        = parts['segment']
        self.observation    = parts['observation']
        self.visit          = parts['visit']

    @property
    def key(self):
        return (self.visit_id, self.exposure)

    @property
    def n_scas(self):
        return len(self.scas)

    @property
    def missing_scas(self):
        return [s for s in range(1, 19) if s not in self.scas]


@dataclass
class DataResults:
    """Everything a downstream tool needs after a MAST list query.

    Carries the raw search results, the unique product list, the
    suffix-filtered subset, and the still-authenticated MastMissions session
    so callers can immediately stream, retrieve, or drill deeper without
    re-authenticating.
    """
    search: dict                    # criteria actually sent to MAST
    data_level: Any                 # 1, 2, 'gw', or None
    results: Any                    # astropy Table from query_criteria
    products: Any                   # unique product list (None if empty)
    filtered: Any                   # suffix-filtered products (or products if data_level=None)
    missions: Any = field(repr=False)   # authenticated MastMissions
    _exposures_cache: Any = field(default=None, repr=False)

    @property
    def n_results(self):
        return len(self.results) if self.results is not None else 0

    @property
    def n_products(self):
        return len(self.filtered) if self.filtered is not None else 0

    @property
    def filenames(self):
        if self.filtered is None or len(self.filtered) == 0:
            return []
        return list(self.filtered['filename'])

    @property
    def exposures(self):
        """List of `Exposure` objects, one per unique (visit_id, exposure_number).

        Built from the metadata `results` table, so each Exposure knows which
        SCAs / filenames belong to it and carries timing metadata for display.
        Mosaic / coadd rows (no per-SCA structure) are skipped.
        """
        if self._exposures_cache is None:
            self._exposures_cache = _group_exposures(self.results, self.data_level)
        return self._exposures_cache

    @property
    def n_exposures(self):
        return len(self.exposures)

    def select(self, key):
        """Return one `Exposure` by index, (visit_id, exp), or Exposure passthrough.

        Accepted keys:
            int              — 1-based index into `res.exposures` (as printed)
            (visit_id, exp)  — tuple of the visit_id string + exposure int
            Exposure         — returned unchanged
        """
        return _select_exposure(self.exposures, key)

    def stream(self, key, *, scas=None, show_progress=True):
        """Stream every SCA of one exposure into memory. See `stream_exposure`."""
        exp = self.select(key)
        return stream_exposure(
            exp, missions=self.missions, scas=scas, show_progress=show_progress,
        )

    def to_fits(self, key, *, out_dir=None, compress=False, sip_degree=4,
                scas=None, overwrite=True, show_progress=True):
        """Stream one exposure and write per-SCA FITS files. Returns the out_dir.

        Convenience wrapper: streams via `stream_exposure`, hands the result
        to `roman_fits.to_fits_files`, then closes the AsdfFiles. Keep the
        streamed buffers around by calling `.stream()` + `to_fits_files()`
        yourself.
        """
        from roman_fits import to_fits_files

        exp = self.select(key)
        af_dict = self.stream(exp, scas=scas, show_progress=show_progress)
        try:
            return to_fits_files(
                af_dict, exp, out_dir=out_dir, compress=compress,
                sip_degree=sip_degree, overwrite=overwrite,
            )
        finally:
            close_streams(af_dict)

    def to_ds9(self, key, *, sip_degree=4, dq_overlay=True,
               with_catalog=True, catalog_radius_arcsec=0.4,
               catalog_color='green', catalog_extended_color='yellow',
               catalog_include_flagged=False, catalog_label_mode='full',
               catalog_show_labels=False,
               ds9_target=None, scas=None, show_progress=True):
        """Stream one exposure and pipe it into DS9 as a WCS mosaic.

        If ``with_catalog=True`` (default), the L4 per-SCA source catalog
        (``cat_sca``) is downloaded for each SCA that has one on MAST and
        drawn as fk5 circle regions on top of the mosaic. Silently skipped
        for SCAs without a catalog.

        Returns the pyds9.DS9 handle so callers can send further XPA commands.
        """
        from roman_fits import to_ds9, download_catalogs

        exp = self.select(key)
        af_dict = self.stream(exp, scas=scas, show_progress=show_progress)
        try:
            catalog_paths = None
            if with_catalog:
                catalog_paths = download_catalogs(
                    exp, self.missions, scas=scas,
                )
            return to_ds9(
                af_dict, exp, sip_degree=sip_degree,
                dq_overlay=dq_overlay,
                catalog_paths=catalog_paths,
                catalog_radius_arcsec=catalog_radius_arcsec,
                catalog_color=catalog_color,
                catalog_extended_color=catalog_extended_color,
                catalog_include_flagged=catalog_include_flagged,
                catalog_label_mode=catalog_label_mode,
                catalog_show_labels=catalog_show_labels,
                ds9_target=ds9_target,
            )
        finally:
            close_streams(af_dict)


def _group_exposures(results, data_level):
    """Group a metadata results Table into a list of Exposure objects."""
    if results is None or len(results) == 0:
        return []

    suffix, ext = DATA_LEVEL_FILE.get(data_level, ('_cal', '.asdf'))

    def _get(row, col, default=None):
        if col not in results.colnames:
            return default
        val = row[col]
        # astropy MaskedColumn: unmasked → return, masked → default
        try:
            if getattr(val, 'mask', False):
                return default
        except Exception:
            pass
        return val

    exposures: dict = {}
    for row in results:
        fsn = _get(row, 'fileSetName')
        if not fsn:
            continue
        m = _FILESET_RE.match(str(fsn))
        if not m:
            # Not a per-SCA row (e.g. mosaic tile) — skip.
            continue

        visit_id = m.group('visit_id')
        exp_num  = int(m.group('exposure'))
        sca      = int(m.group('sca'))
        key = (visit_id, exp_num)

        if key not in exposures:
            exposures[key] = Exposure(
                visit_id=visit_id,
                exposure=exp_num,
                optical_element=_get(row, 'optical_element'),
                exposure_start_time=_get(row, 'exposure_start_time'),
                exposure_time=_get(row, 'exposure_time'),
            )
        exp = exposures[key]
        if sca not in exp.scas:
            exp.scas.append(sca)
            exp.scas.sort()
            exp.filenames.append(f'{fsn}{suffix}{ext}')

    return sorted(exposures.values(), key=lambda e: (e.visit_id, e.exposure))


# ---------------------------------------------------------------------------
# Selecting + streaming an exposure
# ---------------------------------------------------------------------------

def _select_exposure(exposures, key):
    """Resolve `key` into a single Exposure from `exposures`.

    Accepts an int (1-based, matching print_summary numbering), a
    (visit_id, exposure) tuple, or an Exposure passthrough.
    """
    if isinstance(key, Exposure):
        return key

    if isinstance(key, int):
        if key < 1 or key > len(exposures):
            raise IndexError(
                f"Exposure index {key} out of range 1..{len(exposures)}"
            )
        return exposures[key - 1]

    if isinstance(key, tuple) and len(key) == 2:
        visit_id, exp_num = str(key[0]), int(key[1])
        for exp in exposures:
            if exp.visit_id == visit_id and exp.exposure == exp_num:
                return exp
        raise KeyError(
            f"No exposure with visit_id={visit_id!r}, exposure={exp_num}"
        )

    raise TypeError(
        f"Unsupported exposure key {key!r}. Expected int (1-based index), "
        "(visit_id, exposure) tuple, or Exposure."
    )


def stream_exposure(exposure, missions, *, scas=None, show_progress=True):
    """Stream every SCA of `exposure` from MAST into memory.

    Nothing hits the local disk. Uses `MastMissions.read_product` under the
    hood (which is astroquery's authenticated streaming reader).

    Parameters
    ----------
    exposure : Exposure
        The exposure to stream. Use `DataResults.select()` to obtain one, or
        pass `DataResults.stream(key)` and skip this function directly.
    missions : MastMissions
        Authenticated session (typically `res.missions`).
    scas : iterable of int, optional
        Restrict to a subset of SCAs. Default None → every SCA the exposure has.
    show_progress : bool
        Show a tqdm bar over the SCAs.

    Returns
    -------
    dict[int, AsdfFile]
        Mapping SCA integer → open AsdfFile. Call `close_streams()` on the
        result when done, or use `rdm.open(af)` on each entry to get a
        Roman datamodel.
    """
    try:
        from tqdm.auto import tqdm
    except ImportError:  # pragma: no cover
        def tqdm(iterable, **_):    # type: ignore
            return iterable

    if scas is not None:
        wanted = set(int(s) for s in scas)
        pairs = [(s, f) for s, f in zip(exposure.scas, exposure.filenames)
                 if s in wanted]
        missing = wanted - {s for s, _ in pairs}
        if missing:
            _log(f"WARNING: exposure has no filenames for SCAs {sorted(missing)}")
    else:
        pairs = list(zip(exposure.scas, exposure.filenames))

    _log(f"Streaming exposure visit_id={exposure.visit_id} "
         f"exp={exposure.exposure} ({len(pairs)} SCA files)")

    asdf_files = {}
    iterator = tqdm(pairs, desc="Streaming SCAs", disable=not show_progress)
    try:
        for sca, filename in iterator:
            try:
                asdf_files[sca] = missions.read_product(filename)
            except Exception as e:
                _log(f"ERROR streaming SCA {sca:02d} ({filename}): {e}")
                asdf_files[sca] = None
    except KeyboardInterrupt:
        _log(f"Interrupted after {len(asdf_files)} SCAs; closing partial files")
        close_streams(asdf_files)
        raise

    return asdf_files


def close_streams(asdf_files):
    """Close every AsdfFile in a dict returned by `stream_exposure`."""
    for af in asdf_files.values():
        if af is not None and hasattr(af, 'close'):
            try:
                af.close()
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Main entrypoint
# ---------------------------------------------------------------------------

def list_data(
    *,
    program=None,
    execution_plan=None,
    pass_=None,
    segment=None,
    observation=None,
    visit=None,
    detector=None,
    visit_id=None,
    exposure=None,
    optical_element=None,
    exposure_type=None,
    product_type=None,
    sca_only=False,
    data_level=2,
    kinds=None,
    columns=None,
    missions=None,
    token=None,
    enumerate_products=False,
):
    """Query MAST for Roman WFI data matching the given filters.

    All filters are optional; anything left as None matches everything MAST
    returns for that field. See the module docstring for the full reference.

    Parameters
    ----------
    program, pass_, detector, visit_id, exposure, optical_element, exposure_type
        Filter criteria (see module docstring).
    product_type : 'l2' | 'p_visit_coadd' | None
        Server-side filter on which row family MAST returns. Left unset,
        both per-SCA (`l2`) and mosaic-tile (`p_visit_coadd`) rows come back.
    sca_only : bool
        Shortcut for product_type='l2'.
    data_level : 1 | 2 | 'gw' | None
        BACKWARDS-COMPAT shortcut for `kinds`. 1 → ['uncal'], 2 → ['cal'],
        'gw' → ['gw'], None → all kinds. Ignored if `kinds` is set.
    kinds : str | list of str | 'all', optional
        Which product kinds (files per fileSetName row) to synthesize. See
        `PRODUCT_KINDS` for the full menu — 'uncal', 'cal', 'wcs', 'gw',
        'segm_sca', 'cat_sca' (per-SCA rows) plus 'coadd', 'segm', 'cat'
        (mosaic rows). Kinds are family-aware: 'cal' only expands per-SCA
        rows, 'coadd' only expands mosaic rows — so mixing them in one call
        is safe. Overrides `data_level` if both are given.
    columns : list of str, optional
        Columns to request from MAST. Defaults to DEFAULT_COLUMNS.
    missions : MastMissions, optional
        Pre-authenticated session to reuse. If None, connect() is called.
    token : str, optional
        Explicit MAST token override; only consulted when missions is None.
    enumerate_products : bool
        If False (default, FAST), skip the second MAST round-trip and derive
        filenames locally from the `fileSetName` column via Roman's naming
        convention. If True, use MAST's authoritative product-list endpoint
        (get_unique_product_list + filter_products) — this is a second server
        round-trip that can take a very long time on wide searches. In this
        mode the `kinds` list becomes an OR-filter of file_suffix values.

    Returns
    -------
    DataResults
    """
    if missions is None:
        with _timed("Authenticating with MAST"):
            missions = connect(token=token)

    if sca_only:
        if product_type not in (None, 'l2'):
            raise ValueError(
                f"sca_only=True conflicts with product_type={product_type!r}"
            )
        product_type = 'l2'

    search = _build_search(
        program=program, execution_plan=execution_plan, pass_=pass_,
        segment=segment, observation=observation, visit=visit,
        detector=detector, visit_id=visit_id, exposure=exposure,
        optical_element=optical_element, exposure_type=exposure_type,
        product_type=product_type,
    )
    _log(f"Search criteria: {search or '(none — matching everything)'}")
    if not search:
        _log("WARNING: no filters set — this may return a huge result set "
             "and take a very long time.")

    kind_names = _resolve_kinds(kinds, data_level)
    _log(f"Product kinds requested: {kind_names}")

    with _timed("query_criteria (metadata search)"):
        results = missions.query_criteria(
            **search, select_cols=columns or DEFAULT_COLUMNS,
        )
    _log(f"query_criteria returned {len(results)} rows")

    if len(results) == 0:
        return DataResults(
            search=search, data_level=data_level,
            results=results, products=None, filtered=results,
            missions=missions,
        )

    if enumerate_products:
        # Authoritative but slow — one extra MAST round-trip that can time
        # out on wide searches. Kept for the case where you truly need MAST's
        # own product list (e.g. checking for products we don't know about).
        with _timed(f"get_unique_product_list on {len(results)} rows "
                    "(second MAST round-trip; slow on wide searches)"):
            products = missions.get_unique_product_list(results)
        _log(f"get_unique_product_list returned {len(products)} unique products")

        suffixes = sorted({PRODUCT_KINDS[k]['suffix'] for k in kind_names})
        with _timed(f"filter_products(file_suffix={suffixes})"):
            filtered = missions.filter_products(products, file_suffix=suffixes)
        _log(f"filter_products kept {len(filtered)} of {len(products)} products")

        return DataResults(
            search=search, data_level=data_level,
            results=results, products=products, filtered=filtered,
            missions=missions,
        )

    # Fast path: derive filenames locally from fileSetName. No second MAST call.
    _log("Building filenames locally from fileSetName (fast path; "
         "pass enumerate_products=True to hit MAST's product list instead).")
    filtered = _synthesize_products(results, kind_names)
    _log(f"Synthesized {len(filtered)} filenames from {len(results)} rows "
         f"({len(kind_names)} kind(s))")

    return DataResults(
        search=search, data_level=data_level,
        results=results, products=None, filtered=filtered,
        missions=missions,
    )


# ---------------------------------------------------------------------------
# Fast-path filename synthesis
# ---------------------------------------------------------------------------

def _synthesize_products(results, kind_names):
    """Build a fake products Table with one 'filename' column, locally.

    Roman WFI filenames follow the pattern ``{fileSetName}{suffix}{ext}``,
    e.g. ``r0011401057001001001_0001_wfi04_f062_cal.asdf``. `fileSetName` is
    already in the metadata result, so we don't need to ask MAST what
    products exist — we know.

    Kinds are family-aware: a `cal` kind (product_type='l2') is only
    synthesized for rows whose product_type is 'l2', and a `coadd` kind is
    only synthesized for 'p_visit_coadd' rows. This keeps mixed-family
    searches honest — no invented `_cal.asdf` for a mosaic row.
    """
    from astropy.table import Table

    if 'fileSetName' not in results.colnames:
        raise RuntimeError(
            "query_criteria result has no fileSetName column — cannot build "
            "filenames locally. Retry with enumerate_products=True."
        )
    has_pt = 'product_type' in results.colnames

    filenames = []
    for row in results:
        fsn = row['fileSetName']
        if fsn is None or str(fsn) in ('', '--'):
            continue
        row_pt = str(row['product_type']) if has_pt else None

        for k in kind_names:
            spec = PRODUCT_KINDS[k]
            # Only emit this kind for rows of the matching family; if we
            # didn't get product_type back for some reason, be permissive.
            if has_pt and row_pt != spec['product_type']:
                continue
            filenames.append(f"{fsn}{spec['suffix']}{spec['ext']}")

    return Table({'filename': filenames})


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

# Terminal color helpers. Auto-disabled when stdout is not a TTY (piping,
# redirect, capturing output). Respects NO_COLOR (https://no-color.org) and
# the module-level COLOR override set by --color / --no-color.
COLOR = None   # None → auto; True → force on; False → force off

# 256-color codes for the six visit_id chunks (used in header + visit_id string).
# Softer than the base 8-color palette; readable on both light and dark terms.
_CHUNK_COLORS = {
    'program':        '38;5;33',    # blue
    'execution_plan': '38;5;44',    # cyan
    'pass':           '38;5;40',    # green
    'segment':        '38;5;178',   # yellow / amber
    'observation':    '38;5;170',   # magenta
    'visit':          '38;5;203',   # red-orange
}


def _use_color():
    if COLOR is True:
        return True
    if COLOR is False:
        return False
    if os.getenv('NO_COLOR') is not None:
        return False
    return sys.stdout.isatty()


def _ansi(text, *codes):
    """Wrap `text` with ANSI SGR codes, or return it unchanged if color is off."""
    if not codes or not _use_color():
        return text
    return f"\x1b[{';'.join(codes)}m{text}\x1b[0m"


def _colorize_visit_id(vid: str) -> str:
    """Color a 19-digit visit_id by its P/EP/Pass/Seg/Obs/Vis chunks."""
    if len(vid) != 19 or not vid.isdigit():
        return vid
    return (
        _ansi(vid[0:5],   _CHUNK_COLORS['program'])
        + _ansi(vid[5:7],   _CHUNK_COLORS['execution_plan'])
        + _ansi(vid[7:10],  _CHUNK_COLORS['pass'])
        + _ansi(vid[10:13], _CHUNK_COLORS['segment'])
        + _ansi(vid[13:16], _CHUNK_COLORS['observation'])
        + _ansi(vid[16:19], _CHUNK_COLORS['visit'])
    )


def _pad_visible(text: str, visible_text: str, width: int) -> str:
    """Left-pad `text` so its VISIBLE width (ignoring ANSI) is `width`.

    Needed because f-string `<20` counts ANSI escape bytes as characters.
    """
    padding = max(0, width - len(visible_text))
    return text + (' ' * padding)


def _fmt_int(v, w, changed):
    """Right-aligned int column with bold-on-change / dim-when-repeated."""
    if v is None:
        return '-' * w
    s = f"{v:>{w}}"
    return _ansi(s, '1') if changed else _ansi(s, '2')


def _fmt_scas(n_scas, expected=18):
    """Green if SCA count is nominal, yellow if short, else plain."""
    s = f"{n_scas:>5}"
    if n_scas == expected:
        return _ansi(s, '38;5;40')      # green
    if n_scas < expected:
        return _ansi(s, '38;5;220')     # yellow
    return s


def print_summary(res: DataResults, max_rows: int = 50, show_files: bool = False):
    """Print a compact, human-readable summary of a DataResults.

    Exposures (grouped per visit_id + exposure number) are always shown, since
    that's the natural display unit. Pass ``show_files=True`` to also dump the
    flat filename list at the end.

    Terminal color is on when writing to a TTY and off when redirected. Force
    it with ``--color always`` / ``--no-color`` on the CLI, or the ``COLOR``
    module attribute in Python (True / False / None).
    """
    print("=" * 70)
    print("Roman MAST — list_data")
    print("=" * 70)
    print(f"  Criteria      : {res.search or '(none — matching everything)'}")
    print(f"  Data level    : {res.data_level!r}  "
          f"(suffix filter: {DATA_LEVEL_SUFFIX.get(res.data_level, 'none')})")
    print(f"  Result rows   : {res.n_results}")
    print(f"  Products kept : {res.n_products}")
    print(f"  Exposures     : {res.n_exposures}")

    if res.n_products == 0:
        print("\n  (no products match)")
        return

    if res.n_exposures:
        print()
        # Legend: each label is colored the same as its column of digits in the
        # visit_id below. Learn once, then a glance at the visit_id shows the
        # whole (program / execution / pass / …) breakdown.
        legend = ' / '.join([
            _ansi('Program',     _CHUNK_COLORS['program']),
            _ansi('ExecPlan',    _CHUNK_COLORS['execution_plan']),
            _ansi('Pass',        _CHUNK_COLORS['pass']),
            _ansi('Segment',     _CHUNK_COLORS['segment']),
            _ansi('Observation', _CHUNK_COLORS['observation']),
            _ansi('Visit',       _CHUNK_COLORS['visit']),
        ])
        print(f"  visit_id decoded as PPPPPCCAAASSSOOOVVV → {legend}")
        print()

        # Header — chunk-colored labels above matching columns.
        header_bits = [
            _ansi(f"{'#':>3}", '1'),
            _ansi(f"{'Visit ID':<20}", '1'),
            _ansi(f"{'Prog':>5}", '1', _CHUNK_COLORS['program']),
            _ansi(f"{'EP':>2}",   '1', _CHUNK_COLORS['execution_plan']),
            _ansi(f"{'Pass':>4}", '1', _CHUNK_COLORS['pass']),
            _ansi(f"{'Seg':>3}",  '1', _CHUNK_COLORS['segment']),
            _ansi(f"{'Obs':>3}",  '1', _CHUNK_COLORS['observation']),
            _ansi(f"{'Vis':>3}",  '1', _CHUNK_COLORS['visit']),
            _ansi(f"{'Exp':>4}",   '1'),
            _ansi(f"{'Filter':<7}",'1'),
            _ansi(f"{'SCAs':>5}",  '1'),
            _ansi("Start time",    '1'),
        ]
        print("  " + header_bits[0] + "  " + header_bits[1] + " "
              + " ".join(header_bits[2:8]) + "  "
              + header_bits[8] + " " + header_bits[9] + " " + header_bits[10]
              + "  " + header_bits[11])
        print(_ansi("  " + "-" * 100, '2'))

        prev = None
        for i, exp in enumerate(res.exposures):
            if i >= max_rows:
                print(_ansi(f"  ... and {res.n_exposures - max_rows} more exposures", '2'))
                break

            # Change-detection: bold if the field changed vs previous row,
            # dim if it repeats.
            changed = {
                'program':        prev is None or exp.program        != prev.program,
                'execution_plan': prev is None or exp.execution_plan != prev.execution_plan,
                'pass':           prev is None or exp.pass_          != prev.pass_,
                'segment':        prev is None or exp.segment        != prev.segment,
                'observation':    prev is None or exp.observation    != prev.observation,
                'visit':          prev is None or exp.visit          != prev.visit,
                'exposure':       prev is None or exp.exposure       != prev.exposure,
            }

            start = str(exp.exposure_start_time) if exp.exposure_start_time else ''
            filt  = str(exp.optical_element) if exp.optical_element else ''
            vid_colored = _colorize_visit_id(exp.visit_id)

            row = "  " + _ansi(f"{i+1:>3}", '2') + "  "
            row += _pad_visible(vid_colored, exp.visit_id, 20) + " "
            row += _fmt_int(exp.program,        5, changed['program']) + " "
            row += _fmt_int(exp.execution_plan, 2, changed['execution_plan']) + " "
            row += _fmt_int(exp.pass_,          4, changed['pass']) + " "
            row += _fmt_int(exp.segment,        3, changed['segment']) + " "
            row += _fmt_int(exp.observation,    3, changed['observation']) + " "
            row += _fmt_int(exp.visit,          3, changed['visit']) + "  "
            row += _fmt_int(exp.exposure,       4, changed['exposure']) + " "
            row += f"{filt:<7} "
            row += _fmt_scas(exp.n_scas) + "  "
            row += start
            print(row)

            prev = exp

    # --- Products by kind ------------------------------------------------
    # Group the flat filename list into the PRODUCT_KINDS families so a
    # mixed-kind query (--kinds cal,coadd,cat) actually shows what's there.
    if res.n_products:
        by_kind = _classify_products(res.filenames)
        print()
        print(_ansi("  Products by kind:", '1'))
        print(f"    {'kind':<10} {'suffix':<15} {'count':>6}  {'example':<60}")
        print("    " + "-" * 92)
        for kind_name, entries in by_kind.items():
            if not entries:
                continue
            spec = PRODUCT_KINDS.get(kind_name)
            if spec is not None:
                label = _ansi(f"{kind_name:<10}", _CHUNK_COLORS['pass'])
                suffix = f"{spec['suffix']}{spec['ext']}"
            else:
                label = _ansi(f"{kind_name:<10}", '38;5;244')  # unknown → grey
                suffix = '(unrecognized)'
            example = entries[0]
            if len(example) > 60:
                example = example[:57] + '...'
            print(f"    {label} {suffix:<15} {len(entries):>6}  {example}")

    if show_files:
        print("\n  Files:")
        for i, name in enumerate(res.filenames):
            if i >= max_rows:
                print(f"    ... and {res.n_products - max_rows} more")
                break
            print(f"    {name}")


# Coadd/mosaic filenames start with r{PROGRAM}_p_v..., per-SCA filenames with
# r{19-digit visit_id}_. Used to disambiguate kinds with shared suffixes
# (segm / segm_sca and cat / cat_sca both have '_segm' / '_cat').
_COADD_FILENAME_RE = re.compile(r'^r\d{5}_p_v')


def _filename_family(fname):
    """Return 'p_visit_coadd' for a coadd-style filename, else 'l2'."""
    return 'p_visit_coadd' if _COADD_FILENAME_RE.match(fname) else 'l2'


def _classify_products(filenames):
    """Group filenames by which PRODUCT_KINDS entry they match, preserving order.

    Matching is family-aware — the coadd `_segm.asdf` classifies as 'segm'
    (product_type='p_visit_coadd'), while the per-SCA `_segm.asdf` classifies
    as 'segm_sca' (product_type='l2'), even though the suffixes are identical.
    Files that don't match any kind land in an 'other' bucket so surprises
    stay visible.
    """
    # Longest suffix first — e.g. hypothetical '_cal2' shouldn't be shadowed by '_cal'.
    lookup = sorted(
        ((spec['suffix'], spec['ext'], spec['product_type'], name)
         for name, spec in PRODUCT_KINDS.items()),
        key=lambda t: -len(t[0]),
    )
    buckets = {name: [] for name in PRODUCT_KINDS}
    buckets['other'] = []

    for fname in filenames:
        family = _filename_family(fname)
        matched = None
        for suffix, ext, spec_family, name in lookup:
            if spec_family != family:
                continue
            if fname.endswith(f'{suffix}{ext}'):
                matched = name
                break
        buckets[matched if matched else 'other'].append(fname)
    return buckets


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_int_spec(s):
    """Parse '1' / '1-3' / '1,2,5' / '1-3,7' → sorted list of unique ints.

    None or '' returns None (meaning "all", caller's choice).
    """
    if s is None or str(s).strip() in ('', '*', 'all'):
        return None
    out = set()
    for part in str(s).split(','):
        part = part.strip()
        if '-' in part:
            lo, hi = part.split('-', 1)
            out.update(range(int(lo), int(hi) + 1))
        else:
            out.add(int(part))
    return sorted(out)


def add_list_data_args(parser):
    """Attach the standard list_data filter flags to *parser*.

    Shared by every CLI that needs to look up MAST data (roman_mast, roman_fits, …).
    Keeps flag names / help text in one place.
    """
    parser.add_argument('--program',         type=int, default=None,
                        help='APT program ID, e.g. 114')
    parser.add_argument('--execution-plan',  dest='execution_plan',
                        type=int, default=None,
                        help='Execution plan number within the program')
    parser.add_argument('--pass',            dest='pass_', type=int, default=None,
                        metavar='PASS', help='Pass number within the execution plan')
    parser.add_argument('--segment',         type=int, default=None,
                        help='Segment number within the pass')
    parser.add_argument('--observation',     type=int, default=None,
                        help='Observation number within the segment')
    parser.add_argument('--visit',           type=int, default=None,
                        help='Visit number within the observation')
    parser.add_argument('--detector',        default=None,
                        help="Detector — 'WFI04', 'wfi04', or 4")
    parser.add_argument('--visit-id',        default=None,
                        help='Full 19-digit visit ID or wildcard, e.g. 0011401057*')
    parser.add_argument('--exposure',        type=int, default=None,
                        help='Exposure number (last 4 digits of observation_id)')
    parser.add_argument('--optical-element', default=None,
                        help='Optical element, e.g. F062')
    parser.add_argument('--exposure-type',   default=None,
                        help='Exposure type, e.g. WFI_IMAGE')
    parser.add_argument('--product-type',    default=None,
                        help="MAST product_type: 'l2' (raw per-SCA) or "
                             "'p_visit_coadd' (mosaic tile)")
    parser.add_argument('--sca-only',        action='store_true',
                        help='Shortcut for --product-type l2 (raw per-SCA files only)')
    parser.add_argument('--data-level',      default='2',
                        help="Shortcut: 1=uncal, 2=cal, gw=guide-window, "
                             "none=all kinds (default 2). Use --kinds for "
                             "finer control (mosaic coadds, catalogs, ...).")
    parser.add_argument('--kinds',           default=None,
                        help="Product kinds (comma-separated) — overrides "
                             "--data-level if set. 'all' picks every kind. "
                             "See --list-kinds for the menu.")
    parser.add_argument('--list-kinds',      action='store_true',
                        help='Print the PRODUCT_KINDS registry and exit')
    parser.add_argument('--enumerate-products', action='store_true',
                        help='Hit MAST for the authoritative product list instead '
                             'of synthesizing filenames from fileSetName. Slow on '
                             'wide searches; usually unnecessary.')
    parser.add_argument('--quiet',           action='store_true',
                        help='Suppress the [roman_mast] progress diagnostics')
    parser.add_argument('--color',           choices=['auto', 'always', 'never'],
                        default='auto',
                        help='Terminal color: auto (default; on when stdout '
                             'is a TTY), always, never')
    parser.add_argument('--no-color',        dest='color', action='store_const',
                        const='never', help='Alias for --color never')


def list_data_from_args(args):
    """Build a DataResults from a parsed argparse Namespace (add_list_data_args)."""
    if getattr(args, 'quiet', False):
        globals()['VERBOSE'] = False

    color = getattr(args, 'color', 'auto')
    if color == 'always':
        globals()['COLOR'] = True
    elif color == 'never':
        globals()['COLOR'] = False
    # 'auto' leaves COLOR = None → _use_color() checks the TTY.

    kinds = None
    if getattr(args, 'kinds', None):
        raw = str(args.kinds).strip()
        kinds = raw if raw.lower() in ('all', '*') else [
            k.strip() for k in raw.split(',') if k.strip()
        ]

    return list_data(
        program=args.program,
        execution_plan=args.execution_plan,
        pass_=args.pass_,
        segment=args.segment,
        observation=args.observation,
        visit=args.visit,
        detector=args.detector,
        visit_id=args.visit_id,
        exposure=args.exposure,
        optical_element=args.optical_element,
        exposure_type=args.exposure_type,
        product_type=args.product_type,
        sca_only=args.sca_only,
        data_level=_parse_data_level(args.data_level),
        kinds=kinds,
        enumerate_products=args.enumerate_products,
    )


def print_kinds():
    """Dump the PRODUCT_KINDS registry — used by `--list-kinds`."""
    print("Roman product kinds — pass one or more to --kinds (comma-separated):")
    print()
    print(f"  {'name':<10} {'product_type':<15} {'file suffix':<20} description")
    print("  " + "-" * 76)
    for name, spec in PRODUCT_KINDS.items():
        suffix = f"{spec['suffix']}{spec['ext']}"
        print(f"  {name:<10} {spec['product_type']:<15} {suffix:<20} "
              f"{spec['description']}")
    print()
    print("  Notes:")
    print("    * kinds are family-aware — 'cal' only lists per-SCA rows, "
          "'coadd' only mosaic rows.")
    print("    * combine with --product-type to also restrict which rows MAST "
          "returns server-side.")


def _parse_data_level(s):
    if s is None:
        return 2
    s = str(s).strip().lower()
    if s in ('none', 'all', ''):
        return None
    if s == 'gw':
        return 'gw'
    return int(s)


def _cli():
    import argparse

    p = argparse.ArgumentParser(
        description="List available Roman WFI data on MAST.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # The comm_streaming_example.ipynb search
  python roman_mast.py --program 114 --pass 57 --detector WFI04

  # Everything MAST has for one visit
  python roman_mast.py --visit-id 0011401057001001001

  # Level-1 (uncal) files for a program
  python roman_mast.py --program 114 --data-level 1

  # Raw per-SCA exposures only (drop the p_visit_coadd mosaic tiles)
  python roman_mast.py --program 114 --pass 57 --sca-only

  # Show every product kind (cal / uncal / cat / wcs / segm / ...)
  python roman_mast.py --visit-id 0011401057* --data-level none

  # No filters at all — every Roman product on MAST (may be huge)
  python roman_mast.py
""",
    )
    add_list_data_args(p)
    p.add_argument('--max-rows',        type=int, default=50,
                   help='Max exposures/filenames to print (default 50)')
    p.add_argument('--show-files',      action='store_true',
                   help='Also print the flat filename list after the '
                        'per-exposure summary')

    args = p.parse_args()

    if args.list_kinds:
        print_kinds()
        return

    res = list_data_from_args(args)
    print_summary(res, max_rows=args.max_rows, show_files=args.show_files)


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    _cli()
