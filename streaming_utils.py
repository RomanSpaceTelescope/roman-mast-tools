
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




#%%
def group_file_urls(mast_token, visit_id, exp_num='*', data_level=1):
    """Group the URLs for a given visit ID and (optionally) exposure number for a particular data level (1 by default).

    Args:
        mast_token (str): User MAST API Token
        visit_id (str): Visit ID to query for
        exp_num (int, optional): Exposure number to query for. Defaults to all available for the given Visit ID.
        data_level (int, optional): Data level to query for. Defaults to 1.

    Returns:
        dict: Dictionary containing lists of FileGroup objects with the URLs of the queried results, in order of exposure.
    """
    
    file_extensions = {1:'uncal.asdf',
                  2:'cal.asdf',
                  'gw':'gw.asdf'}

    extension = file_extensions[data_level]

    # Query mast with wildcard to find product names for all SCAs
    # Create MastMissions object and assign mission to 'roman'
    missions = MastMissions(mission='roman')

    # Login to search and retrieve Roman data
    missions.login(token=mast_token)

    # Make a list of column names to return in the search results. The results will not be ordered
    # in this order, so we will re-order later. The fileSetName, which we need for retrieval,
    # is not in this list, but will still be returned.
    col_list = ['fileSetName', 'visit_id', 'program', 'execution_plan', 'pass', 'segment', 'visit', 'observation_id', 
                'optical_element', 'exposure_type', 'instrument_name', 'detector', 'productLevel', 
                'product_type', 'exposure_time', 'exposure_start_time', 'exposure_end_time', 'guide_window_id']

    # Since they removed the observation_exposure value for some reason, we have to use
    # a wild card query on the observation ID instead.
    if exp_num is None or exp_num=='':
        exp_num = '*'

    if exp_num != '*':
        exp_num = f'*{exp_num:04}'

    # Create a dictionary of search criteria
    if data_level != 'gw':
        search = {
                'visit_id':visit_id,
                # 'observation_exposure':exp_num,
                'observation_id':exp_num,
                'productLevel':data_level
                }
    else:
        guide_window_id = visit_id
        if exp_num is not None:
            guide_window_id = guide_window_id + str(exp_num)
        search = {
                'visit_id':visit_id,
                'guide_window_id':guide_window_id,
                'productLevel':1
                }

    # Query with column criteria
    # results = missions.query_criteria(
    #     **search,
    #     select_cols=col_list)
    results = missions.query_criteria(
        **search,
        select_cols=col_list)

    # Get the SCA number from the file set names
    scas_to_process = np.unique([get_sca_num(fileset) for fileset in results['fileSetName']])

    # Iterating by result:
    files = {int(scanum):[] for scanum in scas_to_process}
    
    # If we don't have any results, return an empty dictionary.
    if len(results)==0:
        print(f'No results found for {visit_id} for data level {data_level}')
        return {}

    if data_level == 'gw':
        prev_exp_num = int(results[0]['guide_window_id'][-1])
    else:
        prev_exp_num = results[0]['observation_id']


    for i, result in enumerate(results):

        # Get the file set name and the desired file extension
        fileset = result['fileSetName']

        # Get this result's SCA number from the file set name
        scanum = get_sca_num(fileset)

        if data_level == 'gw':
            new_exp_num = int(result['guide_window_id'][-1])
        else:
            # Getting exposure number for this result
            new_exp_num = result['observation_id']
        

        # If we reach the next exposure and there are missing URLs, we need to backfill. 
        if new_exp_num > prev_exp_num:
            files = fill_missing_exps(files)
            prev_exp_num = new_exp_num

        # Make URL using the file set name an dextension
        if data_level == 'gw':
            gw_file_name = fileset.replace(f'{get_exp_num(fileset):04}', str(new_exp_num))
            file_name = f'{gw_file_name}_{extension}'
        else:
            file_name = f'{fileset}_{extension}'

        url = f"https://mast.stsci.edu/search/roman/api/v0.1/retrieve_product?product_name={fileset}%2F{file_name}"

        # Adding url to files dictionary as a file group
        # files[scanum].append(gf.FileGroup(base=url, files=[url]))
        files[scanum].append(url)


    # Backfilling again if there are missing URLs for the last exposure
    files = fill_missing_exps(files)
    
    return files



def stream_to_buffer(url, mast_token=None, show_progress=True, show_buffer_size=True, position=0, cancel_event=None):
    """Stream data from the given URL into the buffer.

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

def stream_file_group_to_buffer(files, exp_num, mast_token=None):
    """Stream all URLs in the given files dictionary and save in the buffer for the given exposure number.

    Args:
        files (dict): {sca: {real_exp_num: url}} as returned by resolve_query().urls
        exp_num (int): Real exposure number (from the filename) to stream.
        mast_token (str, optional): MAST API Authentication Token (if required). Defaults to None.

    Returns:
        dict: {sca: io.BytesIO} buffer for each SCA at the requested exposure.
    """

    urls_to_stream = []   # list of (scanum, url)
    for scanum, exp_map in files.items():

        # Key by real exposure number
        url = exp_map.get(exp_num)

        # If this is not a URL, skip it
        if not isinstance(url, str) or not url.startswith('https'):
            continue

        urls_to_stream.append((scanum, url))


    streamed_buffers = [None] * len(urls_to_stream) # Pre-fill empty list to keep order

    cancel_event = threading.Event()

    with tqdm(total=len(urls_to_stream), desc="Total Files", position=0, leave=True) as overall_bar:
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            try:
                # Submit tasks with their index (i) acting as the tqdm position
                futures = {
                    executor.submit(stream_to_buffer, url, mast_token, show_buffer_size=False, position=i+1, cancel_event=cancel_event): i
                    for i, (_, url) in enumerate(urls_to_stream)
                }

                for future in concurrent.futures.as_completed(futures):
                    original_index = futures[future] # Get the original order index
                    streamed_buffers[original_index] = future.result() # Save in correct spot
                    overall_bar.update(1)

            except KeyboardInterrupt:
                print("\nKeyboardInterrupt detected. Signalling threads to stop streaming...")
                cancel_event.set()

                # Close successfull buffers
                for b in streamed_buffers:
                    if b is not None:
                        b.close()

                # Cancel future threads
                for future in futures.keys():
                    future.cancel()

                raise # Force the ThreadPoolExecutor to shut down quickly

    # Reporting the total buffer size
    total_size = np.sum([x.getbuffer().nbytes for x in streamed_buffers]) / 1024 / 1024

    if total_size > 1000:
        total_size = total_size / 1024
        units = 'GB'
    else:
        units='MB'

    print(f"Final buffer size: {total_size:.2f} {units}")

    return {scanum: buf for (scanum, _), buf in zip(urls_to_stream, streamed_buffers)}



def close_buffer_streams(buffers):
    """Close all buffer streams returned by stream_file_group_to_buffer.

    Args:
        buffers (dict): Dict of {scanum: io.BytesIO} returned by stream_file_group_to_buffer

    Returns:
        None
    """

    for buffer in buffers.values():
        if isinstance(buffer, io.BytesIO):
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


def get_MAST_token(env_filepath='.env'):
    """Get the MAST API authentication token from a .env file.
    The API token must be set in a .env file as:
    MAST_API_TOKEN=<token>

    Args:
        env_filepath (str, optional): Filepath to the .env file. Defaults to local '.env' file.

    Returns:
        str: MAST API authentication token.   
    """

    # Loading environment variables
    dotenv.load_dotenv(env_filepath)

    # Pulling MAST API token from environment variables.
    mast_token = os.environ.get("MAST_API_TOKEN")

    return mast_token

