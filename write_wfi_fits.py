"""Stream Roman WFI SCAs as FITS HDUs and write them to disk or DS9.

write_wfi_fits: Stream all 18 Roman WFI SCAs and write them as compressed FITS files with SIP WCS headers.
  Output: one FITS file per SCA named sca_NN.fits[.fz] in the directory given by output_dir.
  When more than one exposure is selected (via range or wildcard), a sub-folder exp<NN> is
  created inside output_dir for each exposure.
  Uses FITS tile compression (Rice by default) via CompImageHDU.  DS9, Imviz, and astropy
  all read .fz files transparently.  Typical compression ratio ~3-4x for sky-background-limited
  detector images.
  The SIP approximation is accurate to ~0.1 px across the chip, which is sufficient
  for DS9 mosaicking and most analysis tools.

stream_to_ds9: Stream Roman WFI SCAs directly into DS9 without writing to disk.
  Uses pyds9 to load each SCA as a separate frame in DS9, with full WCS headers.
  Allows real-time interactive inspection of streaming data.
"""

import keyring, keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

import warnings
warnings.filterwarnings('ignore')

import os
import re
import io
import argparse
import csv
import numpy as np
import roman_datamodels as rdm
from astropy.io import fits

from streaming_utils import (
    stream_file_group_to_buffer,
    close_buffer_streams,
)
from export_metadata_csv import flatten_metadata
from query_utils import add_query_args, prompt_query_params, resolve_query

try:
    import pyds9
except ImportError:
    pyds9 = None


def stream_to_ds9(visit_id, exp_spec=1, data_level=2, sip_degree=4, verbose=False, sca_spec=None, ds9_target=None):
    """Stream Roman WFI SCAs directly into DS9 without writing to disk.

    Args:
        visit_id (str): Roman visit ID (e.g., '0012401001001002001')
        exp_spec (int | str | None): Exposure selector — single int, range string
            ('1-3'), comma list ('1,2,5'), or None for all (default 1)
        data_level (int): Data level (1=uncal, 2=cal, default 2)
        sip_degree (int): SIP polynomial degree (default 4)
        verbose (bool): Show detailed output (default False)
        sca_spec (str | None): SCA selector — range string ('1-5'), comma list ('1,3,5'),
            or None for all (default None)
        ds9_target (str | None): DS9 target name for pyds9. If None, uses default DS9 instance.

    Returns:
        pyds9.DS9: The DS9 instance with loaded data.

    Raises:
        ImportError: If pyds9 is not installed.
        ValueError: If MAST_API_TOKEN is not available or DS9 is not running.
    """

    if pyds9 is None:
        raise ImportError(
            "pyds9 is required for stream_to_ds9. Install with: pip install pyds9"
        )

    # Get MAST token
    mast_token = get_MAST_token()
    if not mast_token:
        raise ValueError("MAST_API_TOKEN not found in .env file")

    print("="*80)
    print(f"Streaming WFI to DS9")
    print("="*80)
    print(f"Visit ID:           {visit_id}")
    print(f"Exposure spec:      {exp_spec}")
    print(f"SCA spec:           {sca_spec or 'all'}")
    print(f"Data level:         {data_level}")
    print(f"SIP degree:         {sip_degree}")
    print()

    # --- Query -------------------------------------------------------------------
    q = resolve_query(
        visit_id,
        exp_spec=str(exp_spec) if exp_spec is not None else None,
        data_level=data_level,
        sca_spec=sca_spec
    )

    # Build list of unique (visit_id, exp_num) pairs from the summary.
    visit_exp_pairs = []
    for vid in sorted(q.summary.keys()):
        for exp_num in sorted(q.summary[vid]['exposures']):
            visit_exp_pairs.append((vid, exp_num))

    if len(visit_exp_pairs) > 1:
        raise ValueError(
            f"stream_to_ds9 currently supports a single exposure per call. "
            f"Found {len(visit_exp_pairs)} exposures. Use exp_spec to select a single exposure."
        )

    vid, exp_num = visit_exp_pairs[0]

    # Connect to DS9
    try:
        if ds9_target:
            d = pyds9.DS9(target=ds9_target)
        else:
            d = pyds9.DS9()
    except Exception as e:
        raise ValueError(
            f"Failed to connect to DS9: {e}\n"
            f"Make sure DS9 is running: ds9 &"
        )

    print("Streaming data from MAST...")
    buffer_dict = stream_file_group_to_buffer(q.urls, exp_num=exp_num, missions=q.missions)

    print("Building in-memory mosaic MEF...")
    print("-" * 80)

    # Build two parallel MEFs: one for the science data, one for the DQ (bad-
    # pixel) array.  Both use the same SIP WCS per SCA, so DS9 can stitch them
    # as WCS mosaics and overlay the DQ mosaic as a mask on the data frame.
    # `fits mosaicimage wcs` treats every image extension as a tile and
    # stitches by WCS — the equivalent of `ds9 -mosaic *.fits` on the CLI.
    data_hdulist = fits.HDUList([fits.PrimaryHDU()])
    dq_hdulist   = fits.HDUList([fits.PrimaryHDU()])

    # Guide-star window polygons in world coords, one per SCA, drawn as an
    # fk5 region overlay on the mosaic frame.
    gs_region_lines = []

    n_loaded = 0
    n_dq     = 0
    n_gs     = 0
    for scanum in sorted(buffer_dict):
        buf = buffer_dict.get(scanum)
        if buf is None:
            if verbose:
                print(f'  SCA {scanum:02d}: no buffer, skipping')
            continue

        dm   = rdm.open(buf)
        data = dm.data[...].astype(np.float32)

        # --- Build WCS header ---
        wcs  = dm.meta.wcs
        hdr = wcs.to_fits_sip(
            bounding_box=wcs.bounding_box,
            degree=sip_degree,
        )

        # Do NOT stamp SIMPLE/BITPIX/NAXIS* — astropy sets those from the data
        # when the HDU is built; pre-stamping them can conflict with ImageHDU.
        hdr['EXTNAME'] = f'SCA{scanum:02d}'
        hdr['VISITID'] = (vid, 'Roman visit ID')
        hdr['EXPNUM']  = (exp_num, 'Exposure number within visit')
        hdr['SCANUM']  = (scanum, 'SCA number (1-18)')

        # Copy useful metadata fields if present
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

        data_hdulist.append(fits.ImageHDU(data=data, header=hdr.copy(), name=f'SCA{scanum:02d}'))

        # --- Guide-star window region --------------------------------------
        # meta.guide_star.window_[xy]{start,stop} defines a rectangular window
        # in SCA pixel coords (1-based FITS convention).  Project its four
        # corners through the SCA's gwcs to sky and emit an fk5 polygon so
        # the region lands on the correct tile of the WCS mosaic.
        try:
            gs = dm.meta.guide_star
            x0 = int(gs.window_xstart)
            x1 = int(gs.window_xstop)
            y0 = int(gs.window_ystart)
            y1 = int(gs.window_ystop)
            # gwcs pixel_to_world uses 0-based pixel coords; FITS window_* are
            # 1-based, so subtract 1 for the transform.
            xs = np.array([x0, x1, x1, x0], dtype=float) - 1.0
            ys = np.array([y0, y0, y1, y1], dtype=float) - 1.0
            sky = wcs(xs, ys)  # returns (ra, dec) in degrees
            ra_corners, dec_corners = sky[0], sky[1]
            if np.all(np.isfinite(ra_corners)) and np.all(np.isfinite(dec_corners)):
                coords = ','.join(f'{r:.8f},{d:.8f}'
                                  for r, d in zip(ra_corners, dec_corners))
                gs_region_lines.append(
                    f'polygon({coords}) # color=green width=2 text={{SCA{scanum:02d} GS}}'
                )
                n_gs += 1
        except AttributeError:
            pass
        except Exception as e:
            if verbose:
                print(f'  SCA {scanum:02d}: guide-star region skipped ({e})')

        # --- DQ / bad-pixel mask -------------------------------------------
        # Roman ImageModel exposes `dq` as a per-pixel bitmask (uint32); any
        # non-zero bit means the pixel is flagged as bad.  For a DS9 mask
        # overlay we only need the "is bad?" boolean — cast to uint8 so the
        # mosaic MEF stays small.
        try:
            dq = np.asarray(dm.dq[...])
            dq_mask = (dq != 0).astype(np.uint8)
            dq_hdr = hdr.copy()
            dq_hdr['BUNIT']   = 'flag'
            dq_hdr['CONTENT'] = ('DQ_MASK', 'Non-zero = bad pixel (any DQ bit set)')
            dq_hdulist.append(fits.ImageHDU(data=dq_mask, header=dq_hdr, name=f'DQ{scanum:02d}'))
            n_dq += 1
            print(f'  SCA {scanum:02d}: added to mosaic  ({data.shape[0]}x{data.shape[1]} pixels, '
                  f'{int(dq_mask.sum())} bad pixels)')
        except AttributeError:
            print(f'  SCA {scanum:02d}: added to mosaic  ({data.shape[0]}x{data.shape[1]} pixels, '
                  f'no DQ layer)')

        n_loaded += 1

    print("-" * 80)

    if n_loaded == 0:
        close_buffer_streams(buffer_dict)
        raise ValueError("No SCA data was streamed; nothing to load into DS9.")

    # Serialize the data MEF and pipe it to DS9 via XPA.
    data_buf = io.BytesIO()
    data_hdulist.writeto(data_buf)
    data_bytes = data_buf.getvalue()

    print(f"Piping {len(data_bytes)/1e6:.0f} MB data MEF into DS9 as WCS mosaic...")
    d.set('frame delete all')
    d.set('frame new')
    d.set('fits mosaicimage wcs', data_bytes)

    # --- Overlay DQ as a mask -----------------------------------------------
    # DS9's mask feature draws non-zero pixels of a second image on top of the
    # current frame in a chosen color, aligned by WCS.  This is the on-the-wire
    # equivalent of `File → Open As → Mask` from the DS9 GUI, done as another
    # WCS mosaic so the mask tiles line up with the data tiles.
    if n_dq > 0:
        dq_buf = io.BytesIO()
        dq_hdulist.writeto(dq_buf)
        dq_bytes = dq_buf.getvalue()

        print(f"Piping {len(dq_bytes)/1e6:.1f} MB DQ mask MEF into DS9 as overlay...")
        d.set('mask clear')
        d.set('mask color red')
        d.set('mask transparency 50')
        d.set('mask mark nonzero')
        d.set('fits mask mosaicimage wcs', dq_bytes)

    # --- Guide-star region overlay ------------------------------------------
    # Send the fk5 polygons for all SCAs' guide-star windows as a single
    # DS9 region file over XPA.
    if gs_region_lines:
        region_text = "# Region file format: DS9 version 4.1\nfk5\n" + \
                      "\n".join(gs_region_lines) + "\n"
        print(f"Sending {n_gs} guide-star window regions to DS9...")
        d.set('regions delete all')
        d.set('regions', region_text)

    d.set('zoom to fit')

    # Clean up
    close_buffer_streams(buffer_dict)

    extras = []
    if n_dq:
        extras.append(f'DQ mask overlay ({n_dq} SCAs)')
    if n_gs:
        extras.append(f'{n_gs} guide-star regions')
    suffix = f" with {' + '.join(extras)}" if extras else ''
    print(f"\nLoaded {n_loaded} SCAs as a single WCS mosaic in DS9{suffix}")
    print(f"Visit {vid} / Exposure {exp_num}")

    return d


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

    # Collect all metadata rows across all exposures for top-level CSV
    all_metadata_rows = []

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
        buffer_dict = stream_file_group_to_buffer(q.urls, exp_num=exp_num, missions=q.missions)

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

            # Accumulate for top-level CSV
            all_metadata_rows.extend(metadata_rows)

        print(f'\nFiles in: {exp_dir}/')
        print(f'Open in DS9:  ds9 -mosaic {exp_dir}/sca_*.{ext}')

        # --- Clean up -----------------------------------------------------------
        close_buffer_streams(buffer_dict)

    # --- Write aggregated top-level metadata CSV --------------------------------
    if all_metadata_rows and multi_exposure:
        # Collect all unique keys from metadata rows
        all_keys = set()
        for row in all_metadata_rows:
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

        # Generate top-level CSV filename
        top_csv_file = os.path.join(base_dir, f'metadata_all_level{data_level}.csv')

        print(f"\nWriting aggregated metadata CSV: {top_csv_file}")
        print(f"  Rows (all SCAs): {len(all_metadata_rows)}")
        print(f"  Columns (mnemonics): {len(ordered_keys)}")

        with open(top_csv_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=ordered_keys, restval='')
            writer.writeheader()

            for row in all_metadata_rows:
                # Fill missing values
                for key in ordered_keys:
                    if key not in row:
                        row[key] = ''
                writer.writerow(row)

        print(f"✓ Aggregated metadata written to: {top_csv_file}")

    print(f'\nAll done.')
    if multi_exposure:
        pairs_str = ', '.join([f'{vid}_exp{exp:02d}' for vid, exp in visit_exp_pairs])
        print(f'Exposures: {pairs_str}  →  {base_dir}/')

    return base_dir


def main():
    """Interactive or command-line entry point."""
    parser = argparse.ArgumentParser(
        description='Stream Roman WFI SCAs to disk or DS9.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Write to disk (compressed FITS files)
  python write_wfi_fits.py 0012401001001002001 --exp-num 1
  python write_wfi_fits.py 0012401001001002001 --exp-num 2 --data-level 1 --compress
  python write_wfi_fits.py 0012401001001002001 --exp-range 1-3 --output-dir /tmp/fits

  # Stream to DS9 (requires DS9 running and pyds9 installed)
  python write_wfi_fits.py 0012401001001002001 --exp-num 1 --to-ds9
  python write_wfi_fits.py 0012401001001002001 --exp-num 2 --data-level 1 --to-ds9 --sip-degree 3

  # Interactive mode
  python write_wfi_fits.py
        """)

    add_query_args(parser, visit_wildcard=False, exp_mode='flexible', sca_mode='all')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory for FITS files (disk mode only)')
    parser.add_argument('--sip-degree', type=int, default=None,
                        help='SIP polynomial degree (default: 4)')
    parser.add_argument('--compress', action='store_true',
                        help='Use RICE tile compression (.fits.fz) instead of plain FITS (disk mode only)')
    parser.add_argument('--to-ds9', action='store_true',
                        help='Stream directly to DS9 instead of writing to disk')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed output')

    args = parser.parse_args()

    if args.visit_id is None:
        print("Roman WFI Data Streamer")
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

        mode_str = input("Output mode? (disk/ds9, default disk): ").strip().lower()
        to_ds9 = mode_str == 'ds9'

        if to_ds9:
            verbose_str = input("Verbose output? (y/n, default n): ").strip().lower()
            verbose = verbose_str == 'y'
            output_dir = None
            compress = False
        else:
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
        to_ds9 = args.to_ds9
        compress = args.compress
        verbose = args.verbose
        output_dir = args.output_dir

    try:
        if to_ds9:
            stream_to_ds9(
                visit_id=visit_id,
                exp_spec=exp_spec,
                data_level=data_level,
                sip_degree=sip_degree,
                verbose=verbose,
                sca_spec=sca_spec,
            )
        else:
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
