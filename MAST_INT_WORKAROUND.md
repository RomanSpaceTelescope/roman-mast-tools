# Temporary MAST I&T workaround

**Context:** MAST OPS was cleared in preparation for launch. The MAST Integration & Test
(I&T) instance (`mastint.stsci.edu`) retains a copy of selected commissioning datasets
(e.g. program 114 / MRT-7b). This document describes what was changed to point
`roman_mast.py` at the I&T server, and exactly what to revert when MAST OPS is
repopulated.

---

## What was changed in `roman_mast.py`

### 1. New imports (top of file)
```python
from astroquery.mast import MastMissions, Conf   # Conf added
from requests.adapters import HTTPAdapter as _HTTPAdapter  # new
```

### 2. Module-level default (`MAST_SERVER`)
```python
MAST_SERVER = None  # set to 'https://mastint.stsci.edu' to use I&T
```

### 3. `connect()` — three additions
- `server` parameter (falls back to `MAST_SERVER`)
- Sets `Conf.server` and patches `_service_api_connection.MISSIONS_URL` when a
  non-production server is in use
- Mounts `_HttpsUpgradeAdapter` on `missions._auth_obj.session` (see below)

### 4. `_HttpsUpgradeAdapter` class
The I&T auth server (`auth.mastint.stsci.edu`) is only reachable over HTTPS, but
`mastint.stsci.edu/whoami` redirects to `http://auth.mastint.stsci.edu` (port 80),
which times out from outside STScI. The adapter silently rewrites those redirects
to HTTPS before `urllib3` opens the socket.

### 5. `list_data()` — `server` parameter threaded through
Passed to `connect()` when `missions=None`.

### 6. CLI — `--server` flag
`add_list_data_args()` and `list_data_from_args()` expose `--server` so the I&T URL
can be set from the command line without touching source code.

---

## How to use it

**In Python:**
```python
import roman_mast
roman_mast.MAST_SERVER = 'https://mastint.stsci.edu'

# Token must be an I&T token, not a production token.
# Obtain one at: https://auth.mastint.stsci.edu/token?suggested_name=Astroquery&suggested_scope=mast:exclusive_access
res = roman_mast.list_data(program=120, ...)
```

**From the CLI:**
```bash
python roman_mast.py --server https://mastint.stsci.edu --program 120 --list
```

**I&T token:** stored separately from the production token. The notebook
`comm_streaming_example_int.ipynb` loads it from a `.env` file or
`mast_api_token.txt`. Keep the two tokens in separate env vars / files to avoid
accidentally authenticating production with an I&T credential.

---

## How to revert when MAST OPS is live

1. **Reset the default:** `MAST_SERVER = None` (already the default; just don't set it).
2. **Remove `_HttpsUpgradeAdapter`** — delete the class definition and the
   `missions._auth_obj.session.mount(...)` call inside `connect()`.
3. **Remove `from requests.adapters import HTTPAdapter as _HTTPAdapter`**.
4. **Optionally remove `Conf` from the astroquery import** if nothing else uses it.

The `server` parameter on `connect()` / `list_data()` and the `--server` CLI flag
are harmless to keep — they're no-ops when `MAST_SERVER` is `None` and no `--server`
flag is passed.
