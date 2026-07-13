
#%%
import os
import io
import re
import requests
import numpy as np
from tqdm.auto import tqdm
import dotenv
from functools import partial
import concurrent.futures
import threading

# Use null keyring to avoid DBus/SecretService errors in headless environments.
import keyring
import keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

from astroquery.mast.missions import MastMissions
import roman_datamodels as rdm




#%%
def create_missions(mast_token):
    """Create and authenticate a MastMissions session for Roman data.

    Args:
        mast_token (str): User MAST API Token

    Returns:
        MastMissions: Logged-in MastMissions object with mission='roman'
    """
    missions = MastMissions(mission='roman')
    missions.login(token=mast_token)
    return missions


def read_mast_product(missions, filename, show_progress=False):
    """Stream a Roman data product from MAST into an AsdfFile.

    Delegates to missions.read_product() which handles authentication and streaming.
    This requires astroquery built from PR #3593 or later.

    Args:
        missions (MastMissions): Authenticated MastMissions session
        filename (str): Product filename (e.g., 'r0011401057001001001_0001_wfi04_f062_cal.asdf')
        show_progress (bool): Unused (kept for compatibility); astroquery handles progress internally

    Returns:
        AsdfFile: ASDF file object opened from the streamed data, ready for rdm.open()
    """
    return missions.read_product(filename)


#%%
def group_file_urls(mast_token=None, visit_id=None, exp_num='*', data_level=1, missions=None):
    """Query MAST for Roman data products and group filenames by SCA and exposure.

    Args:
        mast_token (str, optional): User MAST API Token. Required if missions is None.
        visit_id (str): Visit ID to query for
        exp_num (str or int, optional): Exposure number wildcard to query for. Defaults to '*' (all).
        data_level (int, optional): Data level to query for (1=uncal, 2=cal). Defaults to 1.
        missions (MastMissions, optional): Pre-authenticated MastMissions session. If provided, mast_token is ignored.

    Returns:
        tuple: (missions, files) where files is dict: {sca: {exp_num: filename}}
    """

    if missions is None:
        if mast_token is None:
            raise ValueError("Either mast_token or missions must be provided")
        missions = create_missions(mast_token)

    file_extensions = {1: 'uncal.asdf', 2: 'cal.asdf', 'gw': 'gw.asdf'}
    extension = file_extensions[data_level]

    col_list = ['fileSetName', 'visit_id', 'program', 'execution_plan', 'pass', 'segment', 'visit', 'observation_id',
                'optical_element', 'exposure_type', 'instrument_name', 'detector', 'productLevel',
                'product_type', 'exposure_time', 'exposure_start_time', 'exposure_end_time', 'guide_window_id']

    if exp_num is None or exp_num == '':
        exp_num = '*'

    if exp_num != '*':
        exp_num = f'*{exp_num:04}'

    if data_level != 'gw':
        search = {'visit_id': visit_id}
        if exp_num != '*':
            search['observation_id'] = exp_num
    else:
        guide_window_id = visit_id
        if exp_num != '*':
            guide_window_id = guide_window_id + str(exp_num)
        search = {
            'visit_id': visit_id,
            'guide_window_id': guide_window_id,
        }

    results = missions.query_criteria(**search, select_cols=col_list)

    if len(results) == 0:
        print(f'No results found for {visit_id} for data level {data_level}')
        return missions, {}

    products = missions.get_unique_product_list(results)
    filtered = missions.filter_products(products, file_suffix=f'_{extension[:-5]}')

    scas_to_process = np.unique([get_sca_num(fname) for fname in filtered['filename']])

    files = {int(scanum): {} for scanum in scas_to_process}

    for row in filtered:
        filename = row['filename']
        scanum = get_sca_num(filename)
        exp_num_real = get_exp_num(filename)

        if scanum is not None and exp_num_real is not None:
            files[scanum][exp_num_real] = filename

    return missions, files



def stream_to_buffer(url, mast_token=None, show_progress=True, show_buffer_size=True, position=0, cancel_event=None):
    """Stream data from the given URL into the buffer.

    DEPRECATED: Use read_mast_product() with an authenticated MastMissions object instead.

    Args:
        url (str): URL to stream from
        mast_token (str, optional): MAST API Authentication Token (if required). Defaults to None.
        show_progress (bool, optional): When True, show a progress bar for the stream. Defaults to True.
        show_buffer_size (bool, optional): When True, print the size of the streamed buffer when done. Defaults to True.
        position (int, optional): Position to use for the tqdm bar if using multiple bars at once. Defaults to 0.

    Returns:
        io.BytesIO: Buffer containing data streamed from the given URL.
    """    

    # If the mast token is provided, use in headers
    if mast_token is not None:
        # Pass your headers to fsspec using storage_options
        storage_options = {
            'headers': {
                'Content-Type': 'application/json',
                'Authorization': f'token {mast_token}'
            }
        }
    else:
        storage_options={}

    with requests.get(url, stream=True, **storage_options) as r:
        r.raise_for_status()

        total_size = int(r.headers.get("Content-Length", 0))

        buffer = io.BytesIO()

        try:
            if show_progress:
                with tqdm(
                    total=total_size,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc="Streaming",
                    dynamic_ncols=True,
                    position=position,
                    leave=False
                ) as progress:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if cancel_event and cancel_event.is_set():
                                raise InterruptedError("Stream cancelled by user.")
                        
                        buffer.write(chunk)
                        progress.update(len(chunk))
            else:
                for chunk in r.iter_content(chunk_size=1024 * 1024):

                    if cancel_event and cancel_event.is_set():
                        raise InterruptedError("Stream cancelled by user.")

                    buffer.write(chunk)
        except (InterruptedError, KeyboardInterrupt):
            # Close the partial buffer
            buffer.close()
            raise  # Re-raise so the executor knows it failed

        buffer.seek(0)

    if show_buffer_size:
        print(f"Final buffer size: {buffer.getbuffer().nbytes / 1024 / 1024:.2f} MB")
        buffer.seek(0)


    return buffer

def stream_file_group_to_buffer(files, exp_num, mast_token=None, missions=None):
    """Stream all files in the given dictionary for a single exposure.

    Args:
        files (dict): {sca: {real_exp_num: filename}} as returned by resolve_query().urls
        exp_num (int): Real exposure number (from the filename) to stream.
        mast_token (str, optional): MAST API Token. Required if missions is None.
        missions (MastMissions, optional): Pre-authenticated MastMissions session. If provided, mast_token is ignored.

    Returns:
        dict: {sca: AsdfFile} for each SCA at the requested exposure.
    """

    if missions is None:
        if mast_token is None:
            raise ValueError("Either mast_token or missions must be provided")
        missions = create_missions(mast_token)

    asdf_files = {}

    files_to_stream = []
    for scanum, exp_map in files.items():
        filename = exp_map.get(exp_num)
        if isinstance(filename, str):
            files_to_stream.append((scanum, filename))

    with tqdm(total=len(files_to_stream), desc="Streaming files", position=0, leave=True) as pbar:
        for scanum, filename in files_to_stream:
            try:
                af = read_mast_product(missions, filename)
                asdf_files[scanum] = af
                pbar.update(1)
            except KeyboardInterrupt:
                print(f"\nKeyboardInterrupt detected. Closing {len(asdf_files)} already-streamed files...")
                for af in asdf_files.values():
                    if af is not None:
                        af.close()
                raise

    return asdf_files



def close_buffer_streams(buffers):
    """Close all ASDF file or buffer streams returned by stream_file_group_to_buffer.

    Args:
        buffers (dict): Dict of {scanum: AsdfFile|BytesIO} returned by stream_file_group_to_buffer

    Returns:
        None
    """

    for buffer in buffers.values():
        if buffer is not None and hasattr(buffer, 'close'):
            buffer.close()

    return


def fill_missing_exps(files):
    """This function checks the given files dictionary for any missing exposure entries and backfills them with null values.
    Note that this will always backfill at the end and assumes that there are no chronological gaps.

    Args:
        files (dict): Dictionary of urls for each SCA. Any list of less length will be backfilled with None.

    Returns:
        _type_: _description_
    """    

    # Getting a list of the number of exposure for each SCA.
    lengths = [len(files[x]) for x in files]

    # Finding the total number of exposures so far
    max_len = np.max(lengths)

    # Backfilling missing exposures with null values
    for sca in files:
        if len(files[sca]) < max_len:
            for i in range(len(files[sca]), max_len):
                files[sca].append(None)

    return files

def get_sca_num(filename):

    # Perform regex search on the file name
    regex_str = r'wfi(\d{2})'
    search_result = re.search(regex_str, filename)

    # If no result was found, return None
    if search_result is None:
        return

    # Return the SCA number as an int
    return int(search_result.group(1))

def get_exp_num(filename):

    # Perform regex search on the file name
    regex_str = r'\d{19}_(\d{4})_'
    search_result = re.search(regex_str, filename)

    # If no result was found, return None
    if search_result is None:
        return

    # Return the exposure number as an integer
    return int(search_result.group(1))


def get_MAST_token(token_filepath='mast_api_token.txt'):
    """Get the MAST API authentication token from a text file.
    The API token must be stored in the file as a single line containing the token.

    Args:
        token_filepath (str, optional): Filepath to the token file. Defaults to 'mast_api_token.txt'.

    Returns:
        str: MAST API authentication token, or None if file not found or empty.
    """

    try:
        with open(token_filepath, 'r') as f:
            mast_token = f.read().strip()
        return mast_token if mast_token else None
    except FileNotFoundError:
        return None

