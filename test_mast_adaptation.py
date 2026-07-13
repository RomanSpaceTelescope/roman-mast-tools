#!/usr/bin/env python
"""Test script to verify astroquery MAST adaptation.

Uses the same example as comm_streaming_example.ipynb:
- Program 114, Pass 57, Detector WFI04
"""

import os
import sys
import keyring
import keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

import warnings
warnings.filterwarnings('ignore')

import roman_datamodels as rdm
from streaming_utils import create_missions, read_mast_product, group_file_urls
from query_utils import resolve_query

def test_create_missions():
    """Test create_missions() function."""
    print("\n" + "="*70)
    print("TEST 1: create_missions()")
    print("="*70)

    from streaming_utils import get_MAST_token
    token = get_MAST_token()
    if not token:
        print("⚠️  Skipping: MAST_API_TOKEN not found")
        return None

    try:
        missions = create_missions(token)
        print(f"✓ Mission: {missions.mission}")
        print(f"✓ Service: {missions.service}")
        print(f"✓ Successfully created authenticated MastMissions object")
        return missions
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_group_file_urls(missions):
    """Test group_file_urls() with astroquery methods."""
    print("\n" + "="*70)
    print("TEST 2: group_file_urls() (using astroquery methods)")
    print("="*70)

    if missions is None:
        print("⚠️  Skipping: no missions object")
        return None, None

    try:
        # Use a very specific search like the notebook example
        # This searches for program 114, pass 57, detector WFI04 (SCA 4)
        visit_id = '0011401057001001001'  # From the notebook: program 114, pass 57

        print(f"Querying for visit_id={visit_id}, data_level=2")
        missions_obj, files = group_file_urls(missions=missions, visit_id=visit_id, exp_num='', data_level=2)

        if not files:
            print("⚠️  No files found")
            return missions_obj, files

        print(f"✓ Found {len(files)} SCAs with data:")
        for sca in sorted(files.keys()):
            exp_map = files[sca]
            print(f"  SCA {sca:02d}: {len(exp_map)} exposures")
            for exp_num in sorted(exp_map.keys())[:3]:  # Show first 3 exps
                filename = exp_map[exp_num]
                print(f"    Exp {exp_num}: {filename}")
            if len(exp_map) > 3:
                print(f"    ... and {len(exp_map) - 3} more")

        # group_file_urls should return {sca: {exp_num: filename}}
        # Verify structure
        for sca in files:
            assert isinstance(files[sca], dict), f"SCA {sca} should map to dict, got {type(files[sca])}"
            for exp_num in files[sca]:
                filename = files[sca][exp_num]
                assert isinstance(filename, str), f"Filename should be str, got {type(filename)}"
                assert filename.endswith('.asdf'), f"Filename should end in .asdf: {filename}"

        print(f"✓ Structure verified: dict[sca][exp_num] = filename")
        return missions_obj, files

    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return missions, None


def test_read_mast_product(missions, files):
    """Test read_mast_product() wrapper."""
    print("\n" + "="*70)
    print("TEST 3: read_mast_product() (streaming single file)")
    print("="*70)

    if missions is None or not files:
        print("⚠️  Skipping: no missions object or files")
        return None

    try:
        # Pick a file to stream (first SCA, first exposure)
        first_sca = sorted(files.keys())[0]
        first_exp = sorted(files[first_sca].keys())[0]
        filename = files[first_sca][first_exp]

        print(f"Streaming SCA {first_sca}, Exp {first_exp}: {filename}")
        af = read_mast_product(missions, filename)

        print(f"✓ Streamed successfully")
        print(f"✓ Type: {type(af).__name__}")

        # Verify we can open it with roman_datamodels
        dm = rdm.open(af)
        print(f"✓ Opened with roman_datamodels.open()")
        print(f"✓ Data shape: {dm.data.shape}")
        print(f"✓ Data dtype: {dm.data.dtype}")

        # Check metadata
        if hasattr(dm, 'meta') and hasattr(dm.meta, 'observation'):
            obs = dm.meta.observation
            print(f"✓ Observation metadata available:")
            for key in ['program', 'pass', 'observation']:
                if hasattr(obs, key):
                    print(f"    {key}: {getattr(obs, key)}")

        af.close()
        print(f"✓ Closed successfully")
        return af

    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_resolve_query():
    """Test resolve_query() integration."""
    print("\n" + "="*70)
    print("TEST 4: resolve_query() (integrated flow)")
    print("="*70)

    from streaming_utils import get_MAST_token
    token = get_MAST_token()
    if not token:
        print("⚠️  Skipping: MAST_API_TOKEN not found")
        return

    try:
        visit_id = '0011401057001001001'
        print(f"Resolving query for visit_id={visit_id}")

        q = resolve_query(visit_id, exp_spec=None, data_level=2, sca_spec=None)

        print(f"\n✓ ResolvedQuery returned:")
        print(f"  visit_id: {q.visit_id}")
        print(f"  exp_nums: {q.exp_nums}")
        print(f"  data_level: {q.data_level}")
        print(f"  scas: {q.scas}")
        print(f"  urls keys: {list(q.urls.keys())}")
        print(f"  missions: {type(q.missions).__name__ if q.missions else None}")

        # Verify missions object is attached
        assert q.missions is not None, "missions should be attached to ResolvedQuery"
        print(f"✓ missions object attached to ResolvedQuery")

        # Verify urls structure
        assert isinstance(q.urls, dict), "urls should be dict"
        for sca in q.urls:
            assert isinstance(q.urls[sca], dict), f"urls[{sca}] should be dict"
            for exp_num in q.urls[sca]:
                filename = q.urls[sca][exp_num]
                assert isinstance(filename, str), f"filename should be str, got {type(filename)}"

        print(f"✓ urls structure verified: dict[sca][exp_num] = filename")

    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    print("\n" + "="*70)
    print("MAST ADAPTATION VERIFICATION TEST")
    print("="*70)
    print("Testing astroquery native MAST retrieval pattern")
    print("Example: Program 114, Pass 57, Detector WFI04 (from comm_streaming_example.ipynb)")

    # Test 1: create_missions
    missions = test_create_missions()

    # Test 2: group_file_urls (with astroquery methods)
    if missions:
        missions_obj, files = test_group_file_urls(missions)

        # Test 3: read_mast_product
        if files:
            test_read_mast_product(missions_obj, files)

    # Test 4: resolve_query (integrated)
    test_resolve_query()

    print("\n" + "="*70)
    print("✅ ALL TESTS COMPLETED")
    print("="*70)
