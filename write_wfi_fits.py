"""Stream all 18 Roman WFI SCAs and write them as compressed FITS files with SIP WCS headers.

Output: one FITS file per SCA named sca_NN.fits[.fz] in the directory given by output_dir.
When more than one exposure is selected (via range or wildcard), a sub-folder exp<NN> is
created inside output_dir for each exposure.
Uses FITS tile compression (Rice by default) via CompImageHDU.  DS9, Imviz, and astropy
all read .fz files transparently.  Typical compression ratio ~3-4x for sky-background-limited
detector images.

The SIP approximation is accurate to ~0.1 px across the chip, which is sufficient
for DS9 mosaicking and most analysis tools.
"""

import keyring, keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

import warnings
warnings.filterwarnings('ignore')

import os
import re
import argparse
import csv
import numpy as np
import roman_datamodels as rdm
from astropy.io import fits

from streaming_utils import (
    get_MAST_token,
    stream_file_group_to_buffer,
    close_buffer_streams,
)
from export_metadata_csv import flatten_metadata
from query_utils import add_query_args, prompt_query_params, resolve_query


def write_wfi_fits(visit_id, exp_spec=1, data_level=2, output_dir=None, sip_degree=4, compress=False, verbose=False, sca_spec=None):
    """Stream Roman WFI SCAs and write as compressed FITS files with SIP WCS headers.

    Args:
        visit_id (str): Roman visit ID (e.g., '0012401001001002001')
            For wildcard patterns, each matched visit is processed separately.
        exp_spec (int | str | None): Exposure selector — single int, range string
            ('1-3'), comma list ('1,2,5'), or None for all (default 1)
        data_level (int): Data level (1=uncal, 2=cal, default 2)
        output_dir (str): Output directory for FITS files. If None, uses default:
            - Single (visit, exp): wfi_fits_{visitid}_exp{expnum}/
            - Multiple (visit, exp) pairs: wfi_fits_{first_visitid}/v{vid}_exp{num}/
        sip_degree (int): SIP polynomial degree (default 4)
        compress (bool): Whether to use tile compression (default False)
        verbose (bool): Show detailed output (default False)
        sca_spec (str | None): SCA selector — range string ('1-5'), comma list ('1,3,5'),
            or None for all (default None)

    Returns:
        str: Path to base output directory
    """

    # Get MAST token
    mast_token = get_MAST_token()
    if not mast_token:
        raise ValueError("MAST_API_TOKEN not found in .env file")

    print("="*80)
    print(f"Writing WFI FITS Files")
    print("="*80)
    print(f"Visit ID:           {visit_id}")
    print(f"Exposure spec:      {exp_spec}")
    print(f"SCA spec:           {sca_spec or 'all'}")
    print(f"Data level:         {data_level}")
    print(f"Output directory:   {output_dir or '(auto)'}")
    print(f"SIP degree:         {sip_degree}")
    print(f"Compression:        {'RICE_1 (.fits.fz)' if compress else 'None (.fits)'}")
    print()

    # --- Query -------------------------------------------------------------------
    q = resolve_query(
        visit_id,
        exp_spec=str(exp_spec) if exp_spec is not None else None,
        data_level=data_level,
        sca_spec=sca_spec
    )

    # Build list of unique (visit_id, exp_num) pairs from the summary.
    # For wildcard queries, multiple visits may have the same exp_num; each is separate.
    visit_exp_pairs = []
    for vid in sorted(q.summary.keys()):
        for exp_num in sorted(q.summary[vid]['exposures']):
            visit_exp_pairs.append((vid, exp_num))

    multi_exposure = len(visit_exp_pairs) > 1

    if output_dir is None:
        if multi_exposure:
            # Use the first visit ID for the base directory name
            first_visit = visit_exp_pairs[0][0]
            base_dir = f'/efs/roman_it_shared/mrizzo/wfi_fits_{first_visit}'
        else:
            vid, exp_num = visit_exp_pairs[0]
            base_dir = f'/efs/roman_it_shared/mrizzo/wfi_fits_{vid}_exp{exp_num:02d}'
    else:
        base_dir = output_dir

    ext = 'fits.fz' if compress else 'fits'

    # --- Loop over (visit_id, exp_num) pairs ---
    for vid, exp_num in visit_exp_pairs:
        if multi_exposure:
            exp_dir = os.path.join(base_dir, f'v{vid}_exp{exp_num:02d}')
        else:
            exp_dir = base_dir

        if multi_exposure:
            print()
            print("="*80)
            print(f"Visit {vid} / Exposure {exp_num}  →  {exp_dir}/")
            print("="*80)

        print("Streaming data from MAST...")
        buffer_dict = stream_file_group_to_buffer(q.urls, exp_num=exp_num, mast_token=mast_token)

        os.makedirs(exp_dir, exist_ok=True)

        # --- Collect metadata for CSV -------------------------------------------
        metadata_rows = []

        # --- Write FITS and collect metadata ------------------------------------
        print("Writing FITS files...")
        print("-" * 80)

        for scanum in sorted(buffer_dict):
            buf = buffer_dict.get(scanum)
            if buf is None:
                if verbose:
                    print(f'  SCA {scanum:02d}: no buffer, skipping')
                continue

            buf.seek(0)
            dm   = rdm.open(buf)
            data = dm.data[...].astype(np.float32)

            # Extract metadata for CSV
            row = {
                'visit': vid,
                'exposure': exp_num,
                'sca': scanum,
                'data_shape': str(dm.data.shape),
                'data_dtype': str(dm.data.dtype),
                'data_min': float(np.nanmin(data)),
                'data_max': float(np.nanmax(data)),
                'data_mean': float(np.nanmean(data)),
                'data_valid_pixels': int(np.isfinite(data).sum()),
                'data_nan_pixels': int(np.isnan(data).sum()),
            }

            # Flatten all metadata
            if hasattr(dm, 'meta'):
                meta_flat = flatten_metadata(dm.meta, prefix='meta')
                row.update(meta_flat)

            metadata_rows.append(row)

            # --- Build WCS header ---
            # gwcs.to_fits_sip() returns (fits.Header, ndarray-of-residuals).
            # The bounding_box is required; Roman gwcs objects carry it.
            wcs  = dm.meta.wcs
            hdr = wcs.to_fits_sip(
                bounding_box=wcs.bounding_box,
                degree=sip_degree,
            )

            # Promote to a primary HDU header and add basic image keywords.
            hdr['SIMPLE']  = True
            hdr['BITPIX']  = -32          # float32
            hdr['NAXIS']   = 2
            hdr['NAXIS1']  = data.shape[1]
            hdr['NAXIS2']  = data.shape[0]
            hdr['EXTNAME'] = f'SCA{scanum:02d}'
            hdr['VISITID'] = (vid, 'Roman visit ID')
            hdr['EXPNUM']  = (exp_num,  'Exposure number within visit')
            hdr['SCANUM']  = (scanum,   'SCA number (1-18)')

            # Copy a handful of useful metadata fields if present.
            for src, dest, comment in [
                ('meta.exposure.start_time',      'DATE-BEG', 'Exposure start (UTC)'),
                ('meta.exposure.end_time',        'DATE-END', 'Exposure end (UTC)'),
                ('meta.exposure.effective_exposure_time', 'EXPTIME', 'Effective exposure time [s]'),
                ('meta.instrument.detector',      'DETECTOR', 'Detector name'),
                ('meta.instrument.optical_element', 'FILTER', 'Optical element / filter'),
            ]:
                try:
                    val = dm
                    for attr in src.split('.'):
                        val = getattr(val, attr)
                    hdr[dest] = (str(val), comment)
                except AttributeError:
                    pass

            if compress:
                # CompImageHDU requires an empty PrimaryHDU; data + WCS go in the compressed extension.
                comp = fits.CompImageHDU(
                    data=data,
                    header=hdr,
                    compression_type='RICE_1',
                    tile_shape=(256, 256),
                )
                hdul     = fits.HDUList([fits.PrimaryHDU(), comp])
                out_path = os.path.join(exp_dir, f'sca_{scanum:02d}.fits.fz')
            else:
                hdul     = fits.HDUList([fits.PrimaryHDU(data=data, header=hdr)])
                out_path = os.path.join(exp_dir, f'sca_{scanum:02d}.fits')

            hdul.writeto(out_path, overwrite=True)
            print(f'  SCA {scanum:02d}: wrote {out_path}  ({data.shape[0]}x{data.shape[1]}, '
                  f'{os.path.getsize(out_path)/1e6:.0f} MB)')

        print("-" * 80)

        # --- Write metadata CSV -------------------------------------------------
        if metadata_rows:
            # Collect all unique keys from metadata rows
            all_keys = set()
            for row in metadata_rows:
                all_keys.update(row.keys())

            all_keys = sorted(all_keys)

            # Reorder keys with user-specified priority order
            priority_keys = [
                'visit', 'exposure', 'sca',
                'data_shape', 'data_dtype', 'data_min', 'data_max', 'data_mean', 'data_nan_pixels', 'data_valid_pixels',
                'meta.exposure.effective_exposure_time',
                'meta.exposure.start_time',
                'meta.exposure.end_time',
                'meta.exposure.frame_time',
                'meta.exposure.hga_move',
                'meta.exposure.ma_table_id',
                'meta.file_date',
                'meta.pointing.target_aperture',
                'meta.exposure.ma_table_name',
                'meta.rcs.active',
                'meta.rcs.bank',
                'meta.rcs.counts',
                'meta.rcs.electronics',
                'meta.rcs.led',
                'meta.source_catalog.tweakreg_catalog_name',
                'meta.statistics.good_pixel_fraction',
                'meta.statistics.image_median',
                'meta.statistics.image_rms',
                'meta.statistics.zodiacal_light',
                'meta.pointing.target_dec',
                'meta.pointing.target_ra',
            ]

            ordered_keys = [k for k in priority_keys if k in all_keys]
            remaining_keys = sorted([k for k in all_keys if k not in ordered_keys])
            ordered_keys.extend(remaining_keys)

            # Generate CSV filename
            csv_file = os.path.join(exp_dir, f'metadata_exp{exp_num:02d}_level{data_level}.csv')

            print(f"\nWriting metadata CSV: {csv_file}")
            print(f"  Rows (SCAs): {len(metadata_rows)}")
            print(f"  Columns (mnemonics): {len(ordered_keys)}")

            with open(csv_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=ordered_keys, restval='')
                writer.writeheader()

                for row in metadata_rows:
                    # Fill missing values
                    for key in ordered_keys:
                        if key not in row:
                            row[key] = ''
                    writer.writerow(row)

            print(f"✓ Metadata written to: {csv_file}")

        print(f'\nFiles in: {exp_dir}/')
        print(f'Open in DS9:  ds9 -mosaicimage wcs {exp_dir}/sca_*.{ext}')

        # --- Clean up -----------------------------------------------------------
        close_buffer_streams(buffer_dict)

    print(f'\nAll done.')
    if multi_exposure:
        pairs_str = ', '.join([f'{vid}_exp{exp:02d}' for vid, exp in visit_exp_pairs])
        print(f'Exposures: {pairs_str}  →  {base_dir}/')

    return base_dir


def main():
    """Interactive or command-line entry point."""
    parser = argparse.ArgumentParser(
        description='Stream Roman WFI SCAs and write as compressed FITS files with SIP WCS headers.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python write_wfi_fits.py 0012401001001002001 --exp-num 1
  python write_wfi_fits.py 0012401001001002001 --exp-num 2 --data-level 1 --compress
  python write_wfi_fits.py 0012401001001002001 --exp-range 1-3 --output-dir /tmp/fits
  python write_wfi_fits.py 0012401001001002001 --output-dir /tmp/fits --sip-degree 3
  python write_wfi_fits.py  (interactive mode)
        """)

    add_query_args(parser, visit_wildcard=False, exp_mode='flexible', sca_mode='all')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory for FITS files')
    parser.add_argument('--sip-degree', type=int, default=None,
                        help='SIP polynomial degree (default: 4)')
    parser.add_argument('--compress', action='store_true',
                        help='Use RICE tile compression (.fits.fz) instead of plain FITS')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed output')

    args = parser.parse_args()

    if args.visit_id is None:
        print("Roman WFI FITS Writer")
        print("=" * 70)
        params = prompt_query_params(
            visit_wildcard=False,
            exp_mode='flexible',
            sca_mode='all',
            defaults={'visit_id': '', 'exp_num': 1},
        )
        visit_id = params['visit_id']
        if not visit_id:
            print("ERROR: Visit ID is required")
            return
        exp_spec = params['exp_spec'] if params['exp_spec'] else 1
        sca_spec = params['sca_spec']
        data_level = params['data_level']

        sip_degree_str = input("Enter SIP degree (default 4): ").strip()
        sip_degree = int(sip_degree_str) if sip_degree_str else 4

        compress_str = input("Use RICE compression? (y/n, default n): ").strip().lower()
        compress = compress_str == 'y'

        verbose_str = input("Verbose output? (y/n, default n): ").strip().lower()
        verbose = verbose_str == 'y'

        output_dir = None
    else:
        visit_id = args.visit_id
        exp_spec = args.exp_range or (args.exp_num if args.exp_num is not None else 1)
        sca_spec = args.scas
        data_level = args.data_level if args.data_level is not None else 2
        sip_degree = args.sip_degree if args.sip_degree is not None else 4
        compress = args.compress
        verbose = args.verbose
        output_dir = args.output_dir

    try:
        write_wfi_fits(
            visit_id=visit_id,
            exp_spec=exp_spec,
            data_level=data_level,
            output_dir=output_dir,
            sip_degree=sip_degree,
            compress=compress,
            verbose=verbose,
            sca_spec=sca_spec,
        )
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
