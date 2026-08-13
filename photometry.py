import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
import matplotlib.pyplot as plt
from dataclasses import dataclass
import os


class SourcePhotometry:
    """
    Detect sources and measure fluxes in background-subtracted images.
    Supports both aperture photometry and PSF fitting.
    Inspired heavily from the online photutils tutorials.
    """

    def __init__(self, fwhm_pix=2.5, detection_sigma=5.0,
                 aperture_radius_fwhm=2.5, annulus_inner_fwhm=6.0,
                 annulus_outer_fwhm=8.0):
        """
        Parameters
        ----------
        fwhm_pix : float
            Expected source FWHM in pixels.
        detection_sigma : float
            Detection threshold in sigma above background.
        aperture_radius_fwhm : float
            Aperture radius as multiple of FWHM.
        annulus_inner_fwhm : float
            Inner annulus radius as multiple of FWHM.
        annulus_outer_fwhm : float
            Outer annulus radius as multiple of FWHM.
        """
        self.fwhm_pix = fwhm_pix
        self.detection_sigma = detection_sigma
        self.aperture_radius = aperture_radius_fwhm * fwhm_pix
        self.annulus_inner = annulus_inner_fwhm * fwhm_pix
        self.annulus_outer = annulus_outer_fwhm * fwhm_pix

    def detect_sources(self, image, background_rms, background_map=None):
        """
        Detect point sources using DAOStarFinder.

        Parameters
        ----------
        image : np.ndarray
            Background-subtracted image.
        background_rms : np.ndarray or float
            Background RMS map or scalar.
        background_map : np.ndarray, optional
            Background map (for computing threshold if needed).

        Returns
        -------
        sources : astropy.table.Table
            Detected source catalog.
        """
        if np.isscalar(background_rms):
            threshold = self.detection_sigma * background_rms
        else:
            threshold = self.detection_sigma * np.median(background_rms)

        finder = DAOStarFinder(
            fwhm=self.fwhm_pix,
            threshold=threshold,
            sharplo=0.2, sharphi=1.5,
            roundlo=-1.0, roundhi=1.0
        )

        sources = finder(image)

        if sources is not None:
            # Sort by flux (brightest first)
            sources.sort('flux', reverse=True)
            print(f"  Detected {len(sources)} sources above "
                  f"{self.detection_sigma}σ threshold")
        else:
            print("  No sources detected!")
            from astropy.table import Table
            sources = Table()

        return sources

    def aperture_photometry_on_frame(self, image, sources, bkg_rms,
                                     local_bkg_subtraction=True):
        """
        Perform aperture photometry on detected sources.

        Parameters
        ----------
        image : np.ndarray
            Science image (background already subtracted globally).
        sources : astropy.table.Table
            Source catalog with 'xcentroid', 'ycentroid' columns.
        bkg_rms : np.ndarray
            Background RMS map for error estimation.
        local_bkg_subtraction : bool
            If True, also subtract local annulus background.

        Returns
        -------
        phot_table : astropy.table.Table
            Photometry results.
        """
        if len(sources) == 0:
            return None

        positions = np.column_stack([sources['xcentroid'],
                                     sources['ycentroid']])

        # Source aperture
        apertures = CircularAperture(positions, r=self.aperture_radius)

        # Background annulus
        annuli = CircularAnnulus(positions,
                                r_in=self.annulus_inner,
                                r_out=self.annulus_outer)

        # Perform photometry
        phot = aperture_photometry(image, [apertures, annuli])

        # Local background from annulus
        annulus_area = annuli.area
        aperture_area = apertures.area

        # Sigma-clipped median in annulus for each source
        local_bkg_per_pix = np.zeros(len(sources))
        for i, pos in enumerate(positions):
            ann_mask = annuli[i].to_mask(method='center')
            ann_data = ann_mask.multiply(image)
            ann_weights = ann_mask.data
            valid = ann_weights > 0
            if np.any(valid):
                _, local_bkg_per_pix[i], _ = sigma_clipped_stats(
                    ann_data[valid] / ann_weights[valid])

        phot['local_bkg_per_pix'] = local_bkg_per_pix
        phot['local_bkg_total'] = local_bkg_per_pix * aperture_area

        if local_bkg_subtraction:
            phot['flux_bkgsub'] = (phot['aperture_sum_0'] -
                                   phot['local_bkg_total'])
        else:
            phot['flux_bkgsub'] = phot['aperture_sum_0']

        # Error estimation (photon noise + background noise)
        if np.isscalar(bkg_rms):
            bkg_noise_per_pix = bkg_rms
        else:
            # Sample RMS at source positions
            bkg_noise_per_pix = np.array([
                bkg_rms[int(np.clip(pos[1], 0, bkg_rms.shape[0]-1)),
                         int(np.clip(pos[0], 0, bkg_rms.shape[1]-1))]
                for pos in positions
            ])

        # Total noise: sqrt(source_counts + npix * bkg_variance)
        source_counts = np.maximum(phot['flux_bkgsub'], 0)
        phot['flux_err'] = np.sqrt(
            source_counts + aperture_area * bkg_noise_per_pix**2
        )

        # Signal-to-noise ratio
        phot['snr'] = phot['flux_bkgsub'] / phot['flux_err']

        return phot

    def psf_photometry_on_frame(self, image, sources, bkg_rms):
        """
        Perform PSF-fitting photometry using an integrated Gaussian PRF.

        Parameters
        ----------
        image : np.ndarray
            Background-subtracted image.
        sources : astropy.table.Table
            Initial source positions.
        bkg_rms : np.ndarray or float
            Background RMS for weighting.

        Returns
        -------
        psf_phot : astropy.table.Table
            PSF photometry results.
        """
        if len(sources) == 0:
            return None

        from photutils.psf import PSFPhotometry, SourceGrouper
        try:
            from photutils.psf import IntegratedGaussianPRF
        except ImportError:
            from photutils.psf import GaussianPRF as IntegratedGaussianPRF

        # Define the PSF model
        sigma_pix = self.fwhm_pix / 2.355
        psf_model = IntegratedGaussianPRF(sigma=sigma_pix)

        # Set up initial guesses from detection
        init_params = {
            'x_0': sources['xcentroid'],
            'y_0': sources['ycentroid'],
            'flux_0': sources['flux']
        }

        # Grouper for blended sources
        grouper = SourceGrouper(min_separation=2.0 * self.fwhm_pix)

        # Fit size
        fit_shape = (int(5 * self.fwhm_pix) | 1,  # ensure odd
                     int(5 * self.fwhm_pix) | 1)

        # Create the fitter
        from astropy.table import Table

        init_table = Table()
        init_table['x_init'] = sources['xcentroid']
        init_table['y_init'] = sources['ycentroid']
        init_table['flux_init'] = sources['flux']

        psfphot = PSFPhotometry(
            psf_model=psf_model,
            fit_shape=fit_shape,
            grouper=grouper,
            aperture_radius=self.aperture_radius
        )

        psf_result = psfphot(image, init_params=init_table)

        return psf_result

    def display_in_ds9(self, image, sources, photometry=None,
                       frame=1, region_file='sources.reg',
                       title='Roman CRNL Background-Subtracted',
                       cmap='gray', scale='zscale',
                       label_sources=True, show_annuli=True,
                       color_by_snr=True, ds9_target='roman',
                       ds9_path=None, start_ds9=True):
        """
        Display image in DS9 with source regions overlaid.
        Uses command-line xpaget/xpaset (bypasses broken pyds9 XPA bindings).
        """
        import os
        import time
        import shutil
        import subprocess
        import tempfile
        import traceback
        from astropy.io import fits as pyfits

        if sources is None or len(sources) == 0:
            print("No sources to display.")
            return None, None

        # Ensure XPA_METHOD is set
        os.environ['XPA_METHOD'] = 'local'

        # =================================================================
        # Step 1: Write region file
        # =================================================================
        region_file = os.path.abspath(region_file)
        self._write_region_file(sources, photometry, region_file,
                                label_sources, show_annuli, color_by_snr)
        cat_file = region_file.replace('.reg', '.cat')
        self._write_ds9_catalog(sources, photometry, cat_file)

        # =================================================================
        # Step 2: Write image to temporary FITS
        # =================================================================
        tmp_fits = os.path.join(tempfile.gettempdir(), 'roman_ds9_display.fits')
        hdu = pyfits.PrimaryHDU(image.astype(np.float32))
        hdu.header['TITLE'] = title
        hdu.writeto(tmp_fits, overwrite=True)
        print(f"  Temp FITS: {tmp_fits}")

        # =================================================================
        # Step 3: Find binaries
        # =================================================================
        xpaset_bin = shutil.which('xpaset')
        xpaget_bin = shutil.which('xpaget')
        xpaaccess_bin = shutil.which('xpaaccess')

        if ds9_path is None:
            ds9_path = shutil.which('ds9')
            if ds9_path is None:
                ds9_path = '/Applications/SAOImageDS9.app/Contents/MacOS/ds9'

        print(f"  ds9: {ds9_path}")
        print(f"  xpaset: {xpaset_bin}")
        print(f"  xpaget: {xpaget_bin}")

        if not xpaset_bin or not xpaget_bin:
            print("  ERROR: xpaset/xpaget not found on PATH")
            return None, region_file

        # =================================================================
        # Step 4: Launch DS9 if needed
        # =================================================================
        def xpa_check(target):
            """Check if a DS9 target is accessible."""
            r = subprocess.run([xpaaccess_bin, target],
                               capture_output=True, text=True,
                               env={**os.environ, 'XPA_METHOD': 'local'})
            return r.stdout.strip() == 'yes'

        def xpa_set(target, cmd, data=None):
            """Send a command to DS9 via xpaset."""
            env = {**os.environ, 'XPA_METHOD': 'local'}
            if data is not None:
                proc = subprocess.run(
                    [xpaset_bin, target, cmd],
                    input=data, capture_output=True, env=env
                )
            else:
                proc = subprocess.run(
                    [xpaset_bin, '-p', target, cmd],
                    capture_output=True, text=True, env=env
                )
            if proc.returncode != 0:
                stderr = proc.stderr if isinstance(proc.stderr, str) else proc.stderr.decode()
                print(f"    xpaset error ({cmd}): {stderr.strip()}")
                return False
            return True

        def xpa_get(target, cmd):
            """Get data from DS9 via xpaget."""
            env = {**os.environ, 'XPA_METHOD': 'local'}
            proc = subprocess.run(
                [xpaget_bin, target, cmd],
                capture_output=True, text=True, env=env
            )
            return proc.stdout.strip()

        if start_ds9 and not xpa_check(ds9_target):
            print(f"  Launching DS9 (target='{ds9_target}')...")
            subprocess.Popen(
                [ds9_path, '-title', ds9_target],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                env={**os.environ, 'XPA_METHOD': 'local'}
            )

            # Wait for registration
            for i in range(30):
                time.sleep(1)
                if xpa_check(ds9_target):
                    print(f"  DS9 registered after {i+1}s")
                    break
            else:
                print("  ERROR: DS9 never registered after 30s")
                print(f"  Manual fallback:")
                print(f"    ds9 {tmp_fits} -regions {region_file}")
                return None, region_file
        else:
            print(f"  DS9 '{ds9_target}' already running")

        # =================================================================
        # Step 5: Send image and commands via xpaset
        # =================================================================
        print(f"[DS9] Sending image and regions...")

        try:
            # Load FITS file
            xpa_set(ds9_target, f'frame {frame}')
            xpa_set(ds9_target, f'fits {tmp_fits}')
            xpa_set(ds9_target, f'scale {scale}')
            xpa_set(ds9_target, f'cmap {cmap}')
            xpa_set(ds9_target, 'zoom to fit')

            # Load regions
            xpa_set(ds9_target, f'regions load {region_file}')
            xpa_set(ds9_target, 'regions showtext yes')

            print(f"[DS9] Done. {len(sources)} sources overlaid.")
            print(f"\n  Files:")
            print(f"    FITS:    {tmp_fits}")
            print(f"    Regions: {region_file}")
            print(f"    Catalog: {cat_file}")
            print(f"\n  Interact from terminal:")
            print(f"    xpaget {ds9_target} regions")
            print(f"    xpaset -p {ds9_target} scale log")
            print(f"    xpaset -p {ds9_target} cmap heat")

        except Exception as e:
            print(f"  ERROR sending data: {type(e).__name__}: {e}")
            traceback.print_exc()
            return None, region_file

        return ds9_target, region_file

    def _write_region_file(self, sources, photometry, region_file,
                           label_sources, show_annuli, color_by_snr):
        """Write DS9 region file with source apertures and annotations."""
        region_lines = []
        region_lines.append('# Region file format: DS9 version 4.1')
        region_lines.append('global color=green dashlist=8 3 width=1 '
                            'font="helvetica 10 normal roman" select=1 '
                            'highlite=1 dash=0 fixed=0 edit=1 move=1 '
                            'delete=1 include=1 source=1')
        region_lines.append('image')

        for i in range(len(sources)):
            # DS9 uses 1-based pixel coordinates
            x = sources['xcentroid'][i] + 1.0
            y = sources['ycentroid'][i] + 1.0

            # Color by SNR
            color = 'green'
            if color_by_snr and photometry is not None and 'snr' in photometry.colnames:
                snr_val = photometry['snr'][i] if i < len(photometry) else 0
                if snr_val >= 50:
                    color = 'green'
                elif snr_val >= 20:
                    color = 'cyan'
                elif snr_val >= 10:
                    color = 'yellow'
                elif snr_val >= 5:
                    color = 'magenta'
                else:
                    color = 'red'

            # Label
            label_text = ""
            if label_sources:
                parts = [f"#{i+1}"]
                if photometry is not None and i < len(photometry):
                    if 'flux_bkgsub' in photometry.colnames:
                        parts.append(f"F={photometry['flux_bkgsub'][i]:.1f}")
                    elif 'flux_weighted' in photometry.colnames:
                        parts.append(f"F={photometry['flux_weighted'][i]:.1f}")
                    if 'snr' in photometry.colnames:
                        parts.append(f"SNR={photometry['snr'][i]:.1f}")
                    if 'local_bkg_per_pix' in photometry.colnames:
                        parts.append(f"bkg={photometry['local_bkg_per_pix'][i]:.2f}")
                else:
                    # Use detection flux from sources table
                    if 'flux' in sources.colnames:
                        parts.append(f"F={sources['flux'][i]:.1f}")
                label_text = " ".join(parts)

            region_lines.append(
                f'circle({x:.3f},{y:.3f},{self.aperture_radius:.2f}) '
                f'# color={color} width=2 text={{{label_text}}}'
            )

            if show_annuli:
                region_lines.append(
                    f'annulus({x:.3f},{y:.3f},'
                    f'{self.annulus_inner:.2f},{self.annulus_outer:.2f}) '
                    f'# color={color} width=1 dash=1'
                )

        with open(region_file, 'w') as f:
            f.write('\n'.join(region_lines))

        print(f"  Region file: {region_file} ({len(sources)} sources)")


    def _write_ds9_catalog(self, sources, photometry, cat_file):
        """
        Write a DS9-compatible catalog file for the Catalog Tool.

        This creates a tab-separated file that DS9 can load via
        Analysis > Catalogs > Load Catalog.

        Parameters
        ----------
        sources : astropy.table.Table
            Detected sources.
        photometry : astropy.table.Table or None
            Photometry results.
        cat_file : str
            Output filename for the catalog.
        """
        lines = []

        # DS9 catalog header
        lines.append('# DS9 Source Catalog')
        lines.append('# Generated by roman_background_source_fit.py')
        lines.append(f'# {len(sources)} sources')
        lines.append('#')

        # Column definitions for DS9 catalog format
        # DS9 needs special markers for coordinate columns
        if photometry is not None and len(photometry) >= len(sources):
            lines.append('# ID\tx_image\ty_image\tflux\tflux_err\tsnr\t'
                         'local_bkg\tsharpness\troundness')
            lines.append('# ---\t-------\t-------\t----\t--------\t---\t'
                         '---------\t---------\t---------')

            for i in range(len(sources)):
                x = sources['xcentroid'][i] + 1.0  # 1-indexed
                y = sources['ycentroid'][i] + 1.0

                flux = (photometry['flux_bkgsub'][i]
                        if 'flux_bkgsub' in photometry.colnames
                        else photometry.get('flux_weighted', [0])[i]
                        if 'flux_weighted' in photometry.colnames else 0)
                flux_err = (photometry['flux_err'][i]
                            if 'flux_err' in photometry.colnames else 0)
                snr = (photometry['snr'][i]
                       if 'snr' in photometry.colnames else 0)
                bkg = (photometry['local_bkg_per_pix'][i]
                       if 'local_bkg_per_pix' in photometry.colnames else 0)
                sharp = (sources['sharpness'][i]
                         if 'sharpness' in sources.colnames else 0)
                rnd = (sources['roundness1'][i]
                       if 'roundness1' in sources.colnames else 0)

                lines.append(
                    f'{i+1}\t{x:.3f}\t{y:.3f}\t{flux:.2f}\t'
                    f'{flux_err:.2f}\t{snr:.1f}\t{bkg:.4f}\t'
                    f'{sharp:.4f}\t{rnd:.4f}'
                )
        else:
            # Minimal catalog with just positions and detection flux
            lines.append('# ID\tx_image\ty_image\tflux\tsharpness\troundness')
            lines.append('# ---\t-------\t-------\t----\t---------\t---------')

            for i in range(len(sources)):
                x = sources['xcentroid'][i] + 1.0
                y = sources['ycentroid'][i] + 1.0
                flux = sources['flux'][i]
                sharp = (sources['sharpness'][i]
                         if 'sharpness' in sources.colnames else 0)
                rnd = (sources['roundness1'][i]
                       if 'roundness1' in sources.colnames else 0)

                lines.append(
                    f'{i+1}\t{x:.3f}\t{y:.3f}\t{flux:.2f}\t'
                    f'{sharp:.4f}\t{rnd:.4f}'
                )

        with open(cat_file, 'w') as f:
            f.write('\n'.join(lines))

        print(f"  DS9 catalog written: {cat_file}")
