"""
This module contains primitives to generate wavelength
calibration solutions from reduced 1-D spectra.
"""
# ------------------------------------------------------------------------------
import multiprocessing
from pathlib import Path
import time
import traceback
import os
import warnings

from astropy.table import Table
from astropy.time import Time, TimeDelta
from astropy import units
from astroquery.simbad import Simbad
import numpy as np
import pandas as pd
import scipy
from scipy.interpolate import LSQUnivariateSpline, interp1d
from scipy.signal import medfilt

import astrodata
from geminidr.core import Spect
from gempy.adlibrary import dataselect
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

from . import parameters_maroonx_spectrum
from .maroonx_fit import maroonx_fit, set_logger, get_logger
from . import maroonx_utils
from .primitives_maroonx_echelle import MAROONXEchelle
from .maroonx_echellespectrum.maroonxspectrum import MXSpectrum
from .maroonx_echellespectrum.wavelengthsolution import WavelengthSolution

# Monkey patch for barycorrpy compatibility with newer astroquery versions
def _patched_get_stellar_data(name=''):
    """
    Fixed version of barycorrpy.utils.get_stellar_data for newer astroquery.
    Based on solution from: https://github.com/shbhuk/barycorrpy/issues/59
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    warning = []
    
    customSimbad = Simbad()
    customSimbad.add_votable_fields('ra', 'dec', 'pmra', 'pmdec', 'plx_value', 'rvz_radvel')
    
    obj = customSimbad.query_object(name)
    if obj is None:
        raise ValueError(f'ERROR: {name} target not found. Check target name or enter RA,Dec,PMRA,PMDec,Plx,RV,Epoch manually\n\n')
    else:        
        warning += [f'{name} queried from SIMBAD.']

    # Check for masked values
    if all([not x for x in [obj.mask[0][i] for i in obj.colnames]])==False:
        warning += ['Masked values present in queried dataset']

    obj = obj.filled(None)
    
    pos = SkyCoord(ra=obj['ra'], dec=obj['dec'], unit=(u.deg, u.deg))
    ra = pos.ra.value[0]
    dec = pos.dec.value[0]
    pmra = obj['pmra'][0]
    pmdec = obj['pmdec'][0]
    plx = obj['plx_value'][0]
    rv = obj['rvz_radvel'][0] * 1000 #SIMBAD output is in km/s. Converting to m/s
    epoch = 2451545.0
    
    star = {'ra':ra,'dec':dec,'pmra':pmra,'pmdec':pmdec,'px':plx,'rv':rv,'epoch':epoch}
    
    # Fill Masked values with None. Again. 
    for i in star:
        if star[i] > 1e10:
            star[i] = None           
       
    warning += [f'Values queried from SIMBAD are {star}']
    return star, warning

# Monkey patch the original function in barycorrpy
import barycorrpy.utils as bcu
import barycorrpy.barycorrpy as bcp
bcu.get_stellar_data = _patched_get_stellar_data   # update source module
bcp.get_stellar_data = _patched_get_stellar_data   # update alias in barycorrpy.barycorrpy

get_BC_vel = bcp.get_BC_vel
exposure_meter_BC_vel = bcp.exposure_meter_BC_vel

# ------------------------------------------------------------------------------

def _get_calibration_wavecal_path():
    """
    Get the path for the calibration flat file.
    Should probably be deprecated when dragons calib is implemented.
    """
    cwd = Path(os.getcwd())
    return cwd # / 'calibrations' / 'processed_wavecal'

def _get_calibration_wavecal(adinputs):
    """
    Match and return calibration etalon file as astrodata object.
    Should probably be deprecated when dragons calib is implemented.
    """
    # Get the calibration etalons
    calib_path = _get_calibration_wavecal_path()

    arm_tag = 'BLUE' if 'BLUE' in adinputs[0].tags else 'RED'
    etalons = dataselect.select_data(
        list(calib_path.glob('*.fits')), tags=[arm_tag, 'ETALON', 'PREPARED'])
    ad_etalons = [astrodata.open(f) for f in etalons]

    adoutputs = []
    for ad in adinputs:
        science_time = ad.ut_datetime()
        # Find the etalon with minimum time difference
        closest_etalon = min(ad_etalons, 
                           key=lambda etalon: abs((etalon.ut_datetime() - science_time).total_seconds()))
        
        adoutputs.append(closest_etalon)
    return adoutputs


class LogExceptions(object):
    """
    Wraps a function, so that a backtrace is written to logger.
    Used to wrap iterative fit, which is called using multiprocessing.
    """

    def __init__(self, f):
        self.f = f

    def __call__(self, *args, **kwargs):
        try:
            return self.f(*args, **kwargs)
        except Exception as e:
            e.original_traceback = traceback.format_tb(e.__traceback__)
            raise e

@parameter_override
class MaroonXSpectrum(MAROONXEchelle, Spect):
    """
    This class contains primitives to reduce MAROON-X 1-D spectra.
    Code in this class takes already produced 1-D reduced spectra
    and utilizes it to generate mappings from pixel to wavelength
    (dynamic wavelength calibratiosn).
    """

    tagset = {"GEMINI", "MAROONX", "SPECT"}

    def __init__(self, adinputs, **kwargs):
        super(MaroonXSpectrum, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'maroonxdr.maroonx.lookups'
        self._param_update(parameters_maroonx_spectrum)

    def staticWavelengthSolution(self, adinputs=None, **params):
        """
        Load static wavelength solution from lookup file.

        This primitive loads pre-computed static wavelength calibration solutions
        from lookup files and attaches them as an extensions.
        The static solutions provide initial wavelength mappings before dynamic
        refinement using etalon measurements.

        Parameters
        ----------
        adinputs : list of AstroData
            Input AstroData objects with 1D box extracted spectra. Must have
            REDUCED_ORDERS_FIBER_* extensions for each fiber.
        fibers : list of int, optional
            Fiber numbers to load wavelength solutions for. If None, all 
            fibers (1-5) are processed. Default is None.
        suffix : str, optional
            Suffix to append to output filenames. Default is empty string.

        Returns
        -------
        list of AstroData
            Modified AstroData objects with static wavelength solutions stored
            as WLS_STATIC_FIBER_* extensions for each requested fiber. 
            Each extension contains a 2D array with wavelength values (nm) 
            indexed by [order, pixel].
        """

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # requested fibers
        fibers = params.get('fibers')

        for ad in adinputs:
            # Load static wavelength solution from the config
            statwavelength_file = maroonx_utils.get_statwavelength_filename(ad)
            log.info(f'Loading static wavelength file: {statwavelength_file}')

            for fiber in range(1, 6):
                # Set up an initial value for all fibers
                setattr(ad[0], f'WLS_STATIC_FIBER_{fiber}', np.zeros((1, 1)))
                    
            for fiber in fibers:
                # Check if the fiber is present in the data
                orders = getattr(ad[0], f'REDUCED_ORDERS_FIBER_{fiber}')
                if orders.size == 1:
                    # Save an empty array for this fiber
                    wls_data = np.zeros((1, 1))
                else:
                    # Load the static wavelength solution for this fiber
                    wls_dict = maroonx_utils.load_statwls_from_fits(
                        statwavelength_file, 
                        ext_name=f'FIBER_{fiber}',
                        orders=orders)
                    wls_data = np.vstack([wls_dict[str(int(order))] for order in orders])

                setattr(ad[0], f'WLS_STATIC_FIBER_{fiber}', wls_data)

            ad.update_filename(suffix=params['suffix'], strip=True)

        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def getPeaksAndPolynomials(self, adinputs=None, **params):
        """
        Extract etalon peak positions and fit line profile parameters.

        TODO: Rewrite docstring to make reference to maroonx_fit module

        Parameters
        ----------
        adinputs : list of AstroData
            Input AstroData objects with 1D box extracted etalon spectra.
            Must have ETALON tag.
        guess_file : str, optional
            Path to file containing initial guess spectrum for peak positions.
            Default is None.
        fibers : list of int, optional
            Fiber numbers to process. Valid values are 1-5. If None,
            automatically detects fibers with ETALON or LFC fiber type.
            Default is None.
        orders : list of int, optional
            Order numbers to process. If None, all orders are processed.
            Default is None.
        degree_sigma : int, optional
            Polynomial degree for Gaussian sigma variation along order.
            Default is 3.
        degree_width : int, optional
            Polynomial degree for box width variation along order.
            Default is 3.
        use_sigma_lr : bool, optional
            If True, use different polynomial coefficients for left and right
            Gaussian wings. Default is False.
        show_plots : bool, optional
            If True, generate diagnostic plots of etalon line fits.
            Automatically disabled if multithreading is True. Default is False.
        plot_path : str, optional
            Directory path for saving diagnostic plots when show_plots is True.
            Default is empty string (current directory).
        multithreading : bool, optional
            If True, use multiprocessing to parallelize fiber/order fitting.
            Disables show_plots option. Default is True.
        iterations : int, optional
            Maximum number of iterative fitting cycles. Default is 10.
        suffix : str, optional
            Suffix to append to output filenames. Default is empty string.

        Returns
        -------
        list of AstroData
            Modified AstroData objects with two new table extensions:

            - PEAKS: Peak parameters including centroid positions, intensities,
              and widths for each detected etalon line.
            - POLY: Polynomial coefficients describing how box width and
              Gaussian sigma vary across each order.

        Notes
        -----
        The iterative fitting algorithm may skip orders with insufficient flux
        (median < 0.005) or fail for orders with excessive cosmic rays or
        artifacts. Failed fits are logged as warnings.

        Known bad pixels are masked during fitting (e.g., pixel 1943 in order
        122, first 400 pixels in truncated red arm order 94 fiber 5).
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        # Extract parameters from params dict
        guess_file = params.get('guess_file')
        fibers = params.get('fibers')
        orders = params.get('orders')
        degree_sigma = params.get('degree_sigma')
        degree_width = params.get('degree_width')
        use_sigma_lr = params.get('use_sigma_lr')
        show_plots = params.get('show_plots', False)
        plot_path = params.get('plot_path', '')
        multithreading = params.get('multithreading')
        iterations = params.get('iterations')

        start_time = time.time()

        if len(adinputs) == 0:
            raise ValueError("No input files")
        elif len(adinputs) >= 1:
            log.fullinfo(f"Extracting etalon lines from {len(adinputs)} files")

        for ad in adinputs:
            log.fullinfo(f"Extracting peaks from {ad.filename}")

            # See which fibers and orders we are extracting
            if fibers is None:
                log.warning("No fibers specified.  Finding all Etalon and LFC fibers and extracting those.")
                fibers_list = []
                for fiber_num, fiber_type in enumerate(ad.fiber_setup(), start=1):
                    if fiber_type in ['Etalon', 'LFC']:
                        fibers_list.append(fiber_num)
                fibers = tuple(fibers_list)

            if orders is None:
                log.warning("No orders specified.  Extracting all orders.")

            # If multithreading is on, disable plotting
            if multithreading and show_plots:
                log.warning("Multithreading is on.  Disabling plotting.")
                show_plots = False

            errors = []
            results = []
            log.fullinfo(f"Extracting fibers {fibers} from {ad.filename}")

            # Set logger for iterative fit
            #maroonx_fit.set_logger(log)
            set_logger(log)

            # Create worker pools for multithreading
            if multithreading:
                pool = multiprocessing.get_context("fork").Pool()
                log.fullinfo(f"Using {pool._processes} processes for extraction")
                for fiber, order, data, guess in maroonx_utils.load_recordings(
                ad, guess_file, fibers, orders):
                    log.fullinfo(f"{ad.filename} - Fitting fiber {fiber} order {order}")

                    # Remove pixels that are known to be bad
                    ############################
                    if order == 122:
                        # Remnant from the old pipeline, we should move this to the BPM at some point
                        data[1943] = np.nan
                        log.warning(f'Removed pixel 1943 in order 122')

                    if fiber == 5 and order == 94 and len(data) > 4000:
                        data[0:399] = 0
                        log.warning(f'Removed first 400 pixels in truncated order 94 of fiber 5')
                    ############################

                    # Define callback functions to save results and errors for multithreading
                    def save_result(x, fiber=fiber, order=order):
                        x.fiber = fiber
                        x.order = order
                        x.recording_time = 0.0  # TODO extract from data set
                        results.append(x)

                    def error_callback(exc, fiber=fiber, order=order):
                        errors.append((fiber, order, str(exc)))
                        log.warning(f"{ad.filename} - Error extracting fiber {fiber} order {order}: {exc}")

                    if not (np.nanmedian(data) > 0.005):
                        log.warning(f'Skipped order {order} for fiber {fiber} for insufficient flux')
                        continue

                    # Asynchronously run iterative fit using the data yielded by the generator function
                    pool.apply_async(
                        maroonx_fit.iterative_fit,
                        kwds=dict(
                            input_spectrum = data,
                            degree_sigma = degree_sigma,
                            degree_width = degree_width,
                            iterations = iterations,
                            guess_spectrum = guess,
                            fiber="{}_{}".format(fiber, order),
                            plot_path = plot_path,
                            use_sigma_lr = use_sigma_lr,
                            show_plots = show_plots),
                        callback=save_result,
                        error_callback=error_callback
                        )

                pool.close()
                pool.join()

            # If multithreading is off, run in serial
            else:
                log.fullinfo("Using single process for extraction")
                for fiber, order, data, guess in maroonx_utils.load_recordings(
                ad, guess_file, fibers, orders):
                    log.fullinfo(f"{ad.filename} - Fitting fiber {fiber} order {order}")

                    # Remove pixels that are known to be bad
                    ############################
                    if order == 122:
                        # Remnant from the old pipeline, we should move this to the BPM at some point
                        data[1943] = np.nan
                        log.warning(f'Removed pixel 1943 in order 122')

                    if fiber == 5 and order == 94 and len(data) > 4000:
                        data[0:399] = 0
                        log.warning(f'Removed first 400 pixels in truncated order 94 of fiber 5')
                    ############################
                    
                    if not (np.nanmedian(data) > 0.005):
                        log.warning(f'Skipped order {order} for fiber {fiber} for insufficient flux')
                        continue

                    # Run iterative fit using the data yielded by the generator function in serial
                    try:
                        output = maroonx_fit.iterative_fit(
                                input_spectrum = data,
                                degree_sigma = degree_sigma,
                                degree_width = degree_width,
                                iterations = iterations,
                                guess_spectrum = guess,
                                fiber="{}_{}".format(fiber, order),
                                plot_path=plot_path,
                                use_sigma_lr = use_sigma_lr,
                                show_plots = show_plots
                                )

                        # Save the results in a manner anologoous to the save_results callback function
                        output.fiber = fiber
                        output.order = order
                        output.recording_time = 0.0  # TODO extract from data set
                        results.append(output)
                    except Exception as exc:
                        errors.append((fiber, order, str(exc)))
                        log.warning(f"{ad.filename} - Error extracting fiber {fiber} order {order}: {exc}")


            # Record the time taken
            end_time = time.time()
            log.fullinfo(f"Finished extracting etalon lines in {end_time - start_time:.2f} seconds")

            if results:
                # Add the results to adinputs
                log.fullinfo(f"Adding peaks extension to {ad.filename}")
                peaks = maroonx_fit.insert_peak_parameters(results)
                log.fullinfo(f"Adding poly extension to {ad.filename}")
                poly = maroonx_fit.insert_polynomial_parameters(results)
                peaks_astropy = Table.from_pandas(peaks)
                poly_astropy = Table.from_pandas(poly)
                ad[0].PEAKS = peaks_astropy
                ad[0].POLY = poly_astropy

            if errors:
                log.warning("Errors were encountered during extraction")
                for fiber, order, msg in errors:
                    log.warning(f"Error extracting fiber {fiber} order {order}: {msg}")

        return adinputs

    def fitAndApplyEtalonWls(self, adinputs=None, **params):
        """
        Calculate spline-based dynamic wavelength solutions from etalon spectra.

        This primitive processes etalon calibration files to create high-precision
        wavelength solutions for instrumental drift correction. The workflow is:

        1. Load etalon peak data from PEAKS and POLY extensions
        2. Apply initial wavelength estimates (ThAr or static)
        3. Guess etalon peak order numbers using known etalon parameters
        4. Calculate instrumental drift from deviation of measured peaks from
           theoretical etalon wavelengths
        5. Fit cubic splines with outlier rejection to create smooth wavelength
           solutions for each fiber and order
        6. Store dynamic wavelength arrays and drift measurements

        Parameters
        ----------
        adinputs : list of AstroData
            Input AstroData objects containing 1D extracted etalon spectra with
            PEAKS and POLY extensions. Must have ETALON tag.
        fibers : list of int, optional
            Fiber numbers containing etalon spectra to process. Valid values
            are 1-5. If None, automatically detected from fiber setup (fibers
            with type 'Etalon'). Default is None.
        symmetric_linefits : bool, optional
            If True, use symmetric line profile fitting for etalon peak
            detection. Default is False.
        n_knots : int, optional
            Number of interior knots for cubic spline interpolation of
            wavelength solutions. Higher values give more flexible fits.
            Default is 30.
        thar : bool, optional
            If True, apply ThAr-based wavelength solution for initial peak
            identification. If False, use static wavelength vectors from
            configuration. Default is False.
        ref_file : str, optional
            Path to reference etalon file for relative drift measurement.
            Not currently implemented. Default is None.
        ref_fiber : int, optional
            Reference fiber number when ref_file is provided. Not currently
            implemented. Default is None.
        suffix : str, optional
            Suffix to append to output filenames. Default is empty string.

        Returns
        -------
        list of AstroData
            Modified AstroData objects with new extensions and header keywords:

            Extensions added:
            - WLS_DYNAMIC_FIBER_* (1-5): 2D arrays of wavelength values (nm)
              indexed by [order, pixel] for each fiber
            - PEAK_DATA: Updated table with wavelength assignments and peak
              order numbers

            Header keywords added:
            - DRIFT_FIBER_* (1-5): Measured instrumental drift in m/s for
              each processed fiber

        Raises
        ------
        NotImplementedError
            If ref_file parameter is provided (reference file functionality
            not yet implemented).

        Notes
        -----
        Instrumental drift is calculated as the velocity offset between
        measured etalon peak wavelengths and theoretical wavelengths predicted
        from the etalon gap size and refractive index model.

        Spline fitting includes 3.5-sigma outlier clipping to reject cosmic
        rays and bad pixels. Edge pixels are extrapolated linearly from the
        nearest reliable knots.

        Fiber 5 drift is reported separately as 'inst_drift' when present.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # get the parameters from the config
        fibers = params.get('fibers')
        symmetric_linefits = params.get('symmetric_linefits')
        n_knots = params.get('n_knots')
        thar = params.get('thar')
        ref_file = params.get('ref_file')

        for ad in adinputs:
            if "ETALON" not in ad.tags:
                log.warning(f"File {ad.filename} is not ETALON. Skipping dynamic wavelength solution fitting.")
                continue

            # Load the etalon spectrum
            mx_spectrum = MXSpectrum(ad, etalon_peaks_symmetric=symmetric_linefits)
            log.fullinfo(f'Processing etalon file: {ad.filename}')
            
            # Determine fibers to process if not provided
            if fibers is None:
                fibers = []
                for fiber_num, fiber_type in enumerate(ad.fiber_setup(), start=1):
                    if fiber_type == 'Etalon':
                        fibers.append(fiber_num)
                log.fullinfo(f'Found etalon fibers: {fibers}')
            
            if ref_file is not None:
                # If a reference file and reference fiber is given, offset the etalon
                # position on all fibers by the offset found in the reference fiber.
                # ref_fiber = params.get('ref_fiber')
                # TODO: Ask Andreas what this is supposed to do because currently
                # we do not know how this works in the old pipeline
                raise NotImplementedError("Reference file not implemented yet")

            # Load reference wavelength solution from the config
            refwavelength_file = maroonx_utils.get_refwavelength_filename(ad)
            log.fullinfo(f'Loading reference wavelength file: {refwavelength_file}')
            
            # If chosen, apply ThAr wls to etalon frame
            if thar == True:    
                ref_wavelength = {
                    1: maroonx_utils.load_refwls_from_fits(refwavelength_file, ext_name='FIBER_2'),
                    2: maroonx_utils.load_refwls_from_fits(refwavelength_file, ext_name='FIBER_2'),
                    3: maroonx_utils.load_refwls_from_fits(refwavelength_file, ext_name='FIBER_3'),
                    4: maroonx_utils.load_refwls_from_fits(refwavelength_file, ext_name='FIBER_4'),
                    5: maroonx_utils.load_refwls_from_fits(refwavelength_file, ext_name='FIBER_4')
                }
                for fiber in fibers:
                    wls_solution = WavelengthSolution(**ref_wavelength[fiber])
                    mx_spectrum.spectra[fiber].apply_wavelength_solution(wls_solution)
            else:
                for fiber in fibers:
                    log.fullinfo(f'Apply static wavelength solution for peak number identification in fiber {fiber}.')
                    mx_spectrum.spectra[fiber].apply_wavelength_vector()

            log.fullinfo(f'Apply etalon parameters from file {refwavelength_file}.')
            parameters = maroonx_utils.load_params_from_fits(refwavelength_file, ext_name='PARAMETERS')
            for fiber in fibers:
                mx_spectrum.spectra[fiber].etalon_parameters = parameters

            # ========================================================= OK line

            # Guess the order #s of the etalon peak positions in the measured spectrum
            wl = None
            m  = None
            mf = None
            o  = None
            x  = None

            inst_drift = None

            drifts = {}

            for fiber in fibers:
                log.info(f'Guess Etalon peak numbers for fiber {fiber}')
                peak_data = mx_spectrum.spectra[fiber].guess_peak_numbers(debug=0)

                if wl is not None:
                    wl = np.concatenate((wl, peak_data["WAVELENGTH_BY_THAR"].values))
                    m  = np.concatenate((m, peak_data["M"].values))
                    mf = np.concatenate((mf, peak_data["M_FRACTION"].values))
                    o  = np.concatenate((o, peak_data["ORDER"].values))
                    x  = np.concatenate((x, peak_data["CENTER"].values))
                else:
                    wl = peak_data["WAVELENGTH_BY_THAR"].values
                    m  = peak_data["M"].values
                    mf = peak_data["M_FRACTION"].values
                    o  = peak_data["ORDER"].values
                    x  = peak_data["CENTER"].values

                # Calculate drift for this fiber
                residuals = _fc2min(parameters, peak_data["M"].values, peak_data["WAVELENGTH_BY_THAR"].values) / peak_data[
                    "WAVELENGTH_BY_THAR"].values * 3e8
                bad = np.where(np.abs(residuals-np.nanmedian(residuals)) > 4.0 * np.nanstd(residuals))
                residuals[bad] = np.nan
                if fiber == 5:
                    inst_drift = np.nanmean(residuals)
                    drifts[fiber] = inst_drift
                else:
                    drift = np.nanmean(residuals)
                    drifts[fiber] = drift
                    #drifts = np.append(drifts, drift)
                
            # ==============================================================
            # This should probably be splitted into a separate primitive

            # Create new wls from etalon peaks based on the fitted etalon gap size and dispersion model.
            wave = {}

            new_peak_data = []

            for fiber in fibers:
                log.info(f'Spline fit for fiber {fiber}')
                wavelengths_all = []
                residuals_all = []
                orders_all = []
                xs_all = []

                # get fiber spectra and orders
                spectra = mx_spectrum.spectra[fiber]
                speactra_peaks = spectra.peak_data
                for o in spectra.orders:
                    center = speactra_peaks["CENTER"][o]

                    if center.values[0] < center.values[10]:
                        x = (center.values)
                        y = (_peak_to_wavelength_spline(speactra_peaks["M"][o],
                                                    spectra.etalon_parameters).values)
                    else:
                        x = (center.values)[::-1]
                        y = (_peak_to_wavelength_spline(speactra_peaks["M"][o],
                                                    spectra.etalon_parameters).values)[::-1]
                    knots = np.linspace(np.min(x) + 1, np.max(x) - 1, n_knots)
                    lsq = scipy.interpolate.LSQUnivariateSpline(x, y, knots, k=3)
                    r = (y - lsq(x))/y * 3e8
                    good = np.where(np.abs(r) < 3.5 * np.nanstd(r))
                    lsq = scipy.interpolate.LSQUnivariateSpline(x[good], y[good], knots, k=3, ext=3)

                    xs_all          = np.append(xs_all,x[good],axis=0)
                    orders_all      = np.append(orders_all,np.ones_like(x[good])*o,axis=0)
                    wavelengths_all = np.append(wavelengths_all,lsq(x[good]),axis=0)
                    residuals_all   = np.append(residuals_all,y[good] - lsq(x[good]),axis=0)

                    xx = np.arange(len(spectra.data.wavelength.iloc[0]))
                    wavelengths = lsq(xx)

                    f = scipy.interpolate.interp1d(x[good][-2:], y[good][-2:],fill_value='extrapolate')
                    wavelengths[wavelengths==np.max(wavelengths)] = f(xx[wavelengths==np.max(wavelengths)])
                    f = scipy.interpolate.interp1d(x[good][:2], y[good][:2], fill_value='extrapolate')
                    wavelengths[wavelengths == np.min(wavelengths)] = f(xx[wavelengths == np.min(wavelengths)])

                    if fiber in wave:
                        wave[fiber].update({str(int(o)): wavelengths})
                    else:
                        wave[fiber] = {str(int(o)): wavelengths}

                new_peak_data.append(speactra_peaks.copy())

                # Concatenate wavelengths arrays for each fiber and save them
                wave_arrays = [wave[fiber][o] for o in wave[fiber].keys()]
                setattr(ad[0], f'WLS_DYNAMIC_FIBER_{fiber}', np.stack(wave_arrays))

            for fiber in {1, 2, 3, 4, 5} - set(fibers):
                # If fiber is not in the list of fibers, save an empty array
                setattr(ad[0], f'WLS_DYNAMIC_FIBER_{fiber}', np.zeros((1, 1)))

            # Reset indices before concatenating to avoid conflicts
            # for peak_df in new_peak_data:
            #     peak_df.reset_index(drop=True, inplace=True)

            # Collect and reformat updated peak data
            new_peak_data = pd.concat(new_peak_data).set_index(['FIBER', 'ORDER', 'M'], drop=False)
            new_peak_data.sort_index(level=[0, 1, 2], inplace=True)
            ad[0].PEAK_DATA = Table.from_pandas(new_peak_data)

            # Save meassured drift in header entries
            for fiber in fibers:
                ad[0].hdr[f'DRIFT_FIBER_{fiber}'] = (round(drifts[fiber], 2), "Drift in m/s")
                log.info(f"Drift for fiber {fiber}: {round(drifts[fiber], 2)} m/s")

            ad.update_filename(suffix=params['suffix'], strip=True)

        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def applyWavelengthSolution(self, adinputs=None, **params):
        """
        Apply drift-corrected wavelength solutions to science spectra.

        This primitive corrects for instrumental wavelength drifts between the
        etalon calibration exposure and science exposure by comparing simultaneous
        etalon measurements in both frames. The workflow is:

        1. Load science spectra and corresponding etalon calibration with
           dynamic wavelength solutions
        2. Compare etalon peak positions in reference fiber between science
           and calibration frames to measure pixel shifts
        3. Fit smooth spline models to pixel shift variations across orders
        4. Apply measured shifts to science fiber etalon peaks in calibration
           frame to correct for drift
        5. Fit spline-based wavelength solutions for corrected science fibers
        6. Calculate instrumental drift from reference fiber etalon peaks

        Parameters
        ----------
        adinputs : list of AstroData
            Input AstroData objects containing 1D extracted science spectra with
            PEAKS and POLY extensions from getPeaksAndPolynomials.
        fibers : list of int, optional
            Science fiber numbers to process. Valid values are typically 2, 3, 4
            for MAROON-X science fibers. Default is [2, 3, 4].
        ref_fiber : int, optional
            Fiber number containing simultaneous etalon spectrum for drift
            measurement. Typically fiber 5 for MAROON-X. Default is 5.
        symmetric_linefits : bool, optional
            If True, use symmetric line profile fitting for etalon peak
            detection. Default is False.
        n_knots : int, optional
            Number of interior knots for cubic spline interpolation of
            wavelength solutions. Default is 30.
        etalon_file : list of AstroData, optional
            Corresponding etalon calibration files with dynamic wavelength
            solutions from fitAndApplyEtalonWls. If None, calibration database
            is queried for matching etalon frames. Default is None.
        suffix : str, optional
            Suffix to append to output filenames. Default is empty string.

        Returns
        -------
        list of AstroData
            Modified AstroData objects with new extensions and header keywords:

            Extensions added:
            - WLS_SIMULTANEOUS_FIBER_* (fibers + ref_fiber): 2D arrays of
              drift-corrected wavelength values (nm) indexed by [order, pixel]

            Header keywords added:
            - INSTRUME_DRIFT: Instrumental drift measured from reference
              fiber in m/s
            - RELATIVE_DRIFT: Relative drift between science and calibration
              etalon frames in m/s

        Raises
        ------
        KeyError
            If etalon peak lists cannot be matched between science and
            calibration frames for a specific fiber/order combination.

        Notes
        -----
        Pixel shifts are calculated by matching etalon peak positions between
        science and calibration frames using nearest-neighbor indexing with
        0.5 pixel tolerance. Outliers beyond 3-sigma from median shift are
        rejected before spline fitting.

        Order 94 in red arm receives special treatment due to truncation,
        with first 600 pixels using median shift from remaining pixels.

        If spline fitting fails for reference fiber (insufficient good peaks),
        the algorithm relaxes outlier threshold to 5-sigma and retries. Complete
        failures are logged as warnings but don't halt processing.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Get parameters from config
        fibers = params.get('fibers', [2, 3, 4])
        ref_fiber = params.get('ref_fiber', 5)
        symmetric_linefits = params.get('symmetric_linefits', False)
        n_knots = params.get('n_knots', 30)
        etalons = params.get('etalon_file')

        if etalons is None:
            etalons = _get_calibration_wavecal(adinputs)

        for science_ad, etalon_ad in zip(*gt.make_lists(adinputs, etalons)):
            log.fullinfo(f"Processing: {science_ad.filename} , {etalon_ad.filename}")
            log.fullinfo(f'Etalon reference fiber: {ref_fiber}')

            # Load the science and etalon spectrum
            science = MXSpectrum(science_ad, etalon_peaks_symmetric=symmetric_linefits)
            etalon = MXSpectrum(etalon_ad, etalon_peaks_symmetric=symmetric_linefits)

            # Load reference wavelength solution from the config
            refwavelength_file = maroonx_utils.get_refwavelength_filename(science_ad)
            parameters = maroonx_utils.load_params_from_fits(refwavelength_file, ext_name='PARAMETERS')
            log.debug(f'Apply etalon parameters from file {refwavelength_file}.')
            for fiber in fibers:
                etalon.spectra[fiber].etalon_parameters = parameters
            etalon.spectra[ref_fiber].etalon_parameters = parameters
            science.spectra[ref_fiber].etalon_parameters = parameters

            # Apply wavelength vector in etalon spectra to the etalon peak positions
            for fiber in fibers:
                etalon.spectra[fiber].apply_wavelength_vector()
            etalon.spectra[ref_fiber].apply_wavelength_vector()


            # Calculate offsets in pixel space between science and etalon frame for reference fiber
            shifts = []
            wavelengths = []
            orders = []
            x_refs = []
            xs = []
            splfits = []

            for o in etalon.spectra[ref_fiber].peak_data["CENTER"].index.levels[0]:
                x_ref = science.spectra[ref_fiber].peak_data["CENTER"][o]
                x = etalon.spectra[ref_fiber].peak_data["CENTER"][o]
                try:
                    shift = x.loc[:] - x_ref.loc[:].reindex(x.loc[:].index, method="nearest", tolerance=0.5)
                    shift[np.abs(shift) > (np.abs(np.nanmedian(shift)) + 3 * np.nanstd(shift))] = np.nan
                    log.fullinfo(f'Removed {np.count_nonzero(np.isnan(shift)):.0f} peaks in offset calculation for fiber {ref_fiber} in order {o}')
                    log.fullinfo(f'Found {np.nanmean(shift):.3f} pix offset for fiber {ref_fiber} in order {o}')
                except KeyError:
                    log.warning(f'Could not reference etalon line lists for fiber {ref_fiber} in order {o}')

                if o == 94 and 'RED' in etalon_ad.tags:
                    shift[0:600] = np.nanmedian(shift[600:])
                mask = np.isnan(shift)
                spl = scipy.interpolate.LSQUnivariateSpline(shift.index[~mask], shift.values[~mask], [1000, 2000, 3000], k=3)

                shifts = np.append(shifts, shift.values * etalon.spectra[ref_fiber].peak_data["DISPERSION_MPS"][o])
                wavelengths = np.append(wavelengths, etalon.spectra[ref_fiber].peak_data["WAVELENGTH_BY_THAR"][o])
                orders = np.append(orders, np.ones_like(shift.values, dtype=int) * o)
                x_refs = np.append(x_refs, x_ref.values)
                xs = np.append(xs, x.values)
                splfits = np.append(splfits, spl(shift.index) * etalon.spectra[ref_fiber].peak_data["DISPERSION_MPS"][o])

                # Apply offsets (smoothed with spline) to etalon frame for science fibers
                for fiber in fibers:
                    etalon.spectra[fiber].peak_data["CENTER"][o] = etalon.spectra[fiber].peak_data["CENTER"][o] - spl(
                        etalon.spectra[fiber].peak_data["CENTER"][o].index)

            # Guess the order #s of the etalon peak positions in the measured spectrum of the science fibers in the etalon spectrum
            for fiber in fibers + [ref_fiber]:
                etalon.spectra[fiber].apply_wavelength_vector()
                etalon.spectra[fiber].guess_peak_numbers(debug=0)
            #etalon.spectra[ref_fiber].apply_wavelength_vector()

            # Guess the order #s of the etalon peak positions in the measured spectrum for reference fiber in science spectrum
            science.spectra[ref_fiber].apply_wavelength_vector()
            peak_data = science.spectra[ref_fiber].guess_peak_numbers(debug=0)
            
            residuals = _fc2min(parameters, peak_data["M"].values, peak_data["WAVELENGTH_BY_THAR"].values) / peak_data["WAVELENGTH_BY_THAR"].values * 3e8
            bad = np.where(np.abs(residuals-np.nanmedian(residuals)) > 4.0 * np.nanstd(residuals))
            residuals[bad] = np.nan
            inst_drift = np.nanmean(residuals)

            # Create new wls from etalon peaks based on the given etalon gap size and dispersion model.
            wave = {}
            for fiber in fibers:

                wavelengths_all = []
                residuals_all = []
                orders_all = []
                xs_all = []

                # get fiber spectra and orders
                spectra = etalon.spectra[fiber]
                speactra_peaks = spectra.peak_data

                for o in spectra.orders:
                    center = speactra_peaks["CENTER"][o]

                    if center.values[0] < center.values[10]:
                        x = (center.values)
                        y = (_peak_to_wavelength_spline(etalon.spectra[fiber].peak_data["M"][o],
                                                    etalon.spectra[fiber].etalon_parameters).values)
                    else:
                        x = (center.values)[::-1]
                        y = (_peak_to_wavelength_spline(etalon.spectra[fiber].peak_data["M"][o],
                                                    etalon.spectra[fiber].etalon_parameters).values)[::-1]
                    knots = np.linspace(min(x) + 1, max(x) - 1, n_knots)
                    mask = np.isnan(x)
                    if np.count_nonzero(mask) > 0:
                        log.warning(f'Found {np.count_nonzero(mask)} NANs in order {o} for fiber {fiber} while building the spline wls')
                    lsq = scipy.interpolate.LSQUnivariateSpline(x[~mask], y[~mask], knots, k=3)
                    r = (y - lsq(x))/y * 3e8
                    good = np.where(np.abs(r-np.nanmedian(r)) < 3.5 * np.nanstd(r))
                    lsq = scipy.interpolate.LSQUnivariateSpline(x[good], y[good], knots, k=3, ext=3)

                    xs_all          = np.append(xs_all,x[good],axis=0)
                    orders_all      = np.append(orders_all,np.ones_like(x[good])*o,axis=0)
                    wavelengths_all = np.append(wavelengths_all,lsq(x[good]),axis=0)
                    residuals_all   = np.append(residuals_all,y[good] - lsq(x[good]),axis=0)

                    xx = np.arange(len(etalon.spectra[fiber].data.wavelength[92]))
                    wavelengths = lsq(xx)

                    f = scipy.interpolate.interp1d(x[good][-2:], y[good][-2:],fill_value='extrapolate')
                    wavelengths[wavelengths==np.max(wavelengths)] = f(xx[wavelengths==np.max(wavelengths)])
                    f = scipy.interpolate.interp1d(x[good][:2], y[good][:2], fill_value='extrapolate')
                    wavelengths[wavelengths == np.min(wavelengths)] = f(xx[wavelengths == np.min(wavelengths)])

                    f = f'fiber_{fiber}'
                    if f in wave:
                        wave[f].update({str(o): wavelengths})
                    else:
                        wave[f] = {str(o): wavelengths}

            # Now for the reference fiber in the science spectrum. This should be refactored at some point
            wavelengths_all = []
            residuals_all = []
            orders_all = []
            xs_all = []

            ref_spectra = science.spectra[ref_fiber]

            for o in ref_spectra.peak_data["CENTER"].index.levels[0]:
                if ref_spectra.peak_data["CENTER"][o].values[0] < ref_spectra.peak_data["CENTER"][o].values[
                    10]:
                    x = (ref_spectra.peak_data["CENTER"][o].values)
                    y = (_peak_to_wavelength_spline(ref_spectra.peak_data["M"][o],
                                                ref_spectra.etalon_parameters).values)
                else:
                    x = (ref_spectra.peak_data["CENTER"][o].values)[::-1]
                    y = (_peak_to_wavelength_spline(ref_spectra.peak_data["M"][o],
                                                ref_spectra.etalon_parameters).values)[::-1]
                knots = np.linspace(min(x) + 1, max(x) - 1, n_knots)
                lsq = scipy.interpolate.LSQUnivariateSpline(x, y, knots, k=3)
                r = (y - lsq(x)) / y * 3e8
                good = np.where(np.abs(r-np.nanmedian(r)) < 3.5 * np.nanstd(r))
                try:
                    lsq = scipy.interpolate.LSQUnivariateSpline(x[good], y[good], knots, k=3, ext=3)
                except:
                    log.warning(f'Spline fit failed for reference fiber, order {o} - try again')
                    good = np.where(np.abs(r-np.nanmedian(r)) < 5 * np.nanstd(r))
                    try:
                        lsq = scipy.interpolate.LSQUnivariateSpline(x[good], y[good], knots, k=3, ext=3)
                        log.fullinfo(f'Spline fit succeeded for reference fiber, order {o} with higher outlier threshold')
                    except:
                        log.warning(f'Spline fit failed again for reference fiber, order {o}')

                xs_all = np.append(xs_all, x[good], axis=0)
                orders_all = np.append(orders_all, np.ones_like(x[good]) * o, axis=0)
                wavelengths_all = np.append(wavelengths_all, lsq(x[good]), axis=0)
                residuals_all = np.append(residuals_all, y[good] - lsq(x[good]), axis=0)

                xx = np.arange(len(etalon.spectra[fiber].data.wavelength[92]))
                wavelengths = lsq(xx)

                f = scipy.interpolate.interp1d(x[good][-2:], y[good][-2:], fill_value='extrapolate')
                wavelengths[wavelengths == np.max(wavelengths)] = f(xx[wavelengths == np.max(wavelengths)])
                f = scipy.interpolate.interp1d(x[good][:2], y[good][:2], fill_value='extrapolate')
                wavelengths[wavelengths == np.min(wavelengths)] = f(xx[wavelengths == np.min(wavelengths)])

                f = f'fiber_{ref_fiber}'
                if f in wave:
                    wave[f].update({str(o): wavelengths})
                else:
                    wave[f] = {str(o): wavelengths}

            # Calculate mean drift in reference fiber 
            rel_drift = np.nanmean(shifts)

            # Concatenate wavelengths arrays for each fiber and save them
            all_fibers = fibers + [ref_fiber]
            for fiber in all_fibers:
                f = f'fiber_{fiber}'
                wave_arrays = [wave[f][o] for o in wave[f].keys()]
                setattr(science_ad[0], f'WLS_SIMULTANEOUS_FIBER_{fiber}', np.stack(wave_arrays))

            for fiber in {1, 2, 3, 4, 5} - set(all_fibers):
                # If fiber is not in the list of fibers, save an empty array
                setattr(science_ad[0], f'WLS_SIMULTANEOUS_FIBER_{fiber}', np.zeros((1, 1)))

            # Save meassured drift in header entries
            if inst_drift is not None:
                science_ad[0].hdr[f'INSTRUME_DRIFT'] = (round(inst_drift, 2), "Drift in m/s")
                science_ad[0].hdr[f'RELATIVE_DRIFT'] = (round(rel_drift, 2), "Drift in m/s")
                log.fullinfo(f'Instrument Drift: {inst_drift:.1f} m/s')
                log.fullinfo(f'Relative drift measured in Fiber {ref_fiber}: {rel_drift:.1f} m/s')    
            log.fullinfo(f'Updated wavelength vector in {science_ad.filename}')

        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def combineFibers(self, adinputs=None, **params):
        """
        Combine multiple fiber spectra into single high-SNR spectrum.

        This primitive combines science fiber spectra (typically fibers 2, 3, 4)
        into a single higher signal-to-noise spectrum using inverse variance
        weighted averaging. The algorithm performs the following steps for each
        echelle order:

        1. Scale all fibers to match the flux level of fiber 3 (reference)
        2. Interpolate fibers 2 and 4 onto fiber 3's wavelength grid
        3. Calculate median intensity across fibers at each pixel
        4. Perform kappa-sigma clipping to reject outliers and cosmic rays
        5. Combine using inverse variance (1/error^2) weights
        6. Interpolate over remaining gaps in combined spectrum

        The resulting combined spectrum has improved SNR while rejecting
        discrepant pixels that may result from cosmic rays or bad pixels
        in individual fibers.

        Parameters
        ----------
        adinputs : list of AstroData
            Input AstroData objects with optimally extracted fiber spectra.
            Must contain OPTIMAL_REDUCED_FIBER_* and OPTIMAL_REDUCED_ERR_*
            extensions for fibers being combined, plus WLS_SIMULTANEOUS_FIBER_*
            wavelength solutions.
        combine_fibers : list of int
            Fiber numbers to combine. For MAROON-X, typically [2, 3, 4] for
            three science fibers. Must contain at least 2 fibers.
        symmetric_linefits : bool, optional
            If True, use wavelength solutions from WLS_SIMULTANEOUS_SYM_FIBER_*
            extensions. If False, use WLS_SIMULTANEOUS_FIBER_* extensions.
            Default is False.
        kappa_sigma : float, optional
            Sigma clipping threshold for outlier rejection. Pixels deviating
            more than kappa_sigma * sqrt(error) from the median are clipped.
            Default is 5.0.
        max_clips : int, optional
            Maximum number of pixels to clip per fiber per order. If exceeded,
            kappa_sigma is automatically increased by 0.5 to prevent
            over-aggressive clipping. Increases continue until kappa_sigma
            reaches 10 or clips < max_clips. Default is 100.
        suffix : str, optional
            Suffix to append to output filenames. Default is empty string.

        Returns
        -------
        list of AstroData
            Modified AstroData objects with new combined fiber extensions:

            - OPTIMAL_REDUCED_FIBER_6 (or _7 if symmetric_linefits=True):
              Combined flux spectrum with same shape as input fibers
            - OPTIMAL_REDUCED_ERR_6 (or _7 if symmetric_linefits=True):
              Combined error spectrum (inverse of summed weights)
            - WLS_SIMULTANEOUS_FIBER_6 or WLS_SIMULTANEOUS_SYM_FIBER_7:
              Wavelength solution copied from fiber 3 (reference fiber)

        Notes
        -----
        Fiber 3 is used as the reference wavelength grid because it typically
        has the best wavelength solution in MAROON-X multi-fiber mode.

        Known bad pixels are additionally masked before combination:
        - Blue arm: pixel 196 in all fibers
        - Red arm: pixels 1793-1794 in all fibers

        If more than 1000 NaN values appear in a combined order, a warning
        is logged but processing continues. These typically result from edge
        effects or extended bad pixel regions.

        Flux scaling is relative, so combined spectra preserve relative
        spectral shapes but not absolute flux calibration.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        
        combine_fibers = params["combine_fibers"]
        symmetric_linefits = params["symmetric_linefits"]
        kappa_sigma = params["kappa_sigma"]
        max_clips = params["max_clips"]

        for ad in adinputs:
            
            # Get fiber data
            fiber_data = {}
            fiber_errors = {}
            fiber_wls = {}
            for fib in combine_fibers:
                fiber_data[fib] = getattr(ad[0], f"OPTIMAL_REDUCED_FIBER_{fib}")
                fiber_errors[fib] = getattr(ad[0], f"OPTIMAL_REDUCED_ERR_{fib}")

                if symmetric_linefits:
                    fiber_wls[fib] = getattr(ad[0], f"WLS_SIMULTANEOUS_FIBER_{fib}")
                else:
                    fiber_wls[fib] = getattr(ad[0], f"WLS_SIMULTANEOUS_FIBER_{fib}")
            
            # Get orders from one of the fibers, not very relevant which one
            orders = getattr(ad[0], f"REDUCED_ORDERS_FIBER_{fib}")
            
            combined_intensity = {}
            combined_error = {}
            combined_wavelength = {}
            
            orig_kappa_sigma = kappa_sigma
            
            for order_idx, order in enumerate(orders):
                kappa_sigma = orig_kappa_sigma
                
                intensity1 = fiber_data[2][order_idx]
                intensity2 = fiber_data[3][order_idx]
                scale1 = np.nansum(intensity1 / np.nansum(intensity2))
                intensity1 = intensity1 / scale1
                intensity3 = fiber_data[4][order_idx]
                scale3 = np.nansum(intensity3 / np.nansum(intensity2))
                intensity3 = intensity3 / scale3

                error1 = fiber_errors[2][order_idx]
                error1 = error1 / scale1
                error2 = fiber_errors[3][order_idx]
                error3 = fiber_errors[4][order_idx]
                error3 = error3 / scale3

                wave1 = fiber_wls[2][order_idx]
                wave2 = fiber_wls[3][order_idx]
                wave3 = fiber_wls[4][order_idx]
                
                # Apply instrument-specific masking
                if 'BLUE' in ad.tags:
                    log.fullinfo('Masking additional pixels in blue arm at pixel 196')
                    intensity1[196] = np.nan
                    intensity2[196] = np.nan
                    intensity3[196] = np.nan
                
                if 'RED' in ad.tags:
                    log.fullinfo('Masking additional pixels in red arm at pixels 1793-1794')
                    intensity1[1793:1794] = np.nan
                    intensity2[1793:1794] = np.nan
                    intensity3[1793:1794] = np.nan

                mask1 = np.isnan(intensity1)
                mask2 = np.isnan(intensity2)
                mask3 = np.isnan(intensity3)

                try:
                    m = np.zeros_like(wave1, dtype=bool)
                    m[np.unique(wave1, return_index=True, return_inverse=True)[1]] = True
                    mask1[~m] = True

                    f1 = interp1d(wave1[~mask1], intensity1[~mask1], fill_value='extrapolate', kind='slinear')
                    f1e = interp1d(wave1[~mask1], error1[~mask1], fill_value='extrapolate', kind='slinear')
                    intensity1_2 = f1(wave2)
                    error1_2 = np.abs(f1e(wave2))
                except:
                    log.error(f'Interpolation of flux in science fiber 1 of order {order} failed')
                    intensity1_2 = intensity1
                    error1_2 = np.abs(error1)*np.nan

                try:
                    m = np.zeros_like(wave3, dtype=bool)
                    m[np.unique(wave3, return_index=True, return_inverse=True)[1]] = True
                    mask3[~m] = True

                    f3 = interp1d(wave3[~mask3], intensity3[~mask3], fill_value='extrapolate', kind='slinear')
                    f3e = interp1d(wave3[~mask3], error3[~mask3], fill_value='extrapolate', kind='slinear')
                    intensity3_2 = f3(wave2)
                    error3_2 = np.abs(f3e(wave2))
                except:
                    log.error(f'Interpolation of flux in science fiber 3 of order {order} failed')
                    intensity3_2 = intensity3
                    error3_2 = np.abs(error3)*np.nan


                intensity1_2[mask1] = np.nan
                intensity3_2[mask3] = np.nan

                median_intensity = np.nanmedian([intensity1_2, intensity2, intensity3_2], axis=0)

                weights1_2 = 1.0 / error1_2
                weights2   = 1.0 / error2
                weights3_2 = 1.0 / error3_2

                weights1_2[np.isnan(weights1_2)] = 0
                weights2[np.isnan(weights2)] = 0
                weights3_2[np.isnan(weights3_2)] = 0

                weights1_2[mask1] = 0
                weights2[mask2] = 0
                weights3_2[mask3] = 0

                weights1_2[np.isnan(intensity1_2)] = 0
                weights2[np.isnan(intensity2)] = 0
                weights3_2[np.isnan(intensity3_2)] = 0

                clip1 = np.where(np.nan_to_num(np.abs(intensity1_2 - median_intensity)) > kappa_sigma * np.nan_to_num(np.sqrt(error1_2)))
                clip2 = np.where(np.nan_to_num(np.abs(intensity2 - median_intensity)) > kappa_sigma * np.nan_to_num(np.sqrt(error2)))
                clip3 = np.where(np.nan_to_num(np.abs(intensity3_2 - median_intensity)) > kappa_sigma * np.nan_to_num(np.sqrt(error3_2)))

                clip1_n = np.size(clip1)
                clip2_n = np.size(clip2)
                clip3_n = np.size(clip3)

                while(max([clip1_n,clip2_n,clip3_n]) > max_clips and kappa_sigma<10):
                    log.warning(f'Number of maximum clipped pixels exceeded in order {order}')
                    kappa_sigma = kappa_sigma + 0.5

                    clip1 = np.where(
                        np.nan_to_num(np.abs(intensity1_2 - median_intensity)) > kappa_sigma * np.nan_to_num(np.sqrt(error1_2)))
                    clip2 = np.where(
                        np.nan_to_num(np.abs(intensity2 - median_intensity)) > kappa_sigma * np.nan_to_num(np.sqrt(error2)))
                    clip3 = np.where(
                        np.nan_to_num(np.abs(intensity3_2 - median_intensity)) > kappa_sigma * np.nan_to_num(np.sqrt(error3_2)))

                    clip1_n = np.size(clip1)
                    clip2_n = np.size(clip2)
                    clip3_n = np.size(clip3)

                log.info(f'Kappa_sigma in order {order}: {kappa_sigma:.1f}')
                log.info(f'Clipped {clip1_n} pixels in fiber 2 of order {order}')
                log.info(f'Clipped {clip2_n} pixels in fiber 3 of order {order}')
                log.info(f'Clipped {clip3_n} pixels in fiber 4 of order {order}')

                weights1_2[clip1] = 0
                weights2[clip2] = 0
                weights3_2[clip3] = 0

                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore', r'All-NaN slice encountered')

                    weights = np.nansum([weights1_2, weights2, weights3_2], axis=0)
                    bad = np.where(weights == 0)
                    weights1_2[bad] = np.nan
                    weights2[bad] = np.nan
                    weights3_2[bad] = np.nan

                    intensity = np.average([
                        np.nan_to_num(intensity1_2), 
                        np.nan_to_num(intensity2), 
                        np.nan_to_num(intensity3_2)],
                        weights=[weights1_2, weights2, weights3_2],
                        axis=0)
                    error = np.abs(1.0 / np.sum([weights1_2, weights2, weights3_2], axis=0))

                mask = np.isnan(intensity)
                if 1000 > np.count_nonzero(mask) >= 0:
                    log.info(f'Number of NANs fixed in order {order}: {np.count_nonzero(mask)}')
                    f = interp1d(wave2[~mask], intensity[~mask], fill_value='extrapolate', kind='slinear')
                    intensity = f(wave2)
                else:
                    log.warning(f'Too many NANs found in order {order}: {np.count_nonzero(mask)}')

                error[mask] = 1e6
                error[np.isnan(error)] = 1e6

                combined_wavelength.update({order: wave2})
                combined_intensity.update({order: intensity})
                combined_error.update({order: error})


            # ===================================================================================

            # Store combined results as new extensions
            combined_intensity_array = np.vstack([combined_intensity[o] for o in orders])
            combined_error_array = np.vstack([combined_error[o] for o in orders])
            combined_wavelength_array = np.vstack([combined_wavelength[o] for o in orders])

            # Define final combined fiber number
            target_fiber = 6 if not symmetric_linefits else 7

            setattr(ad[0], f"OPTIMAL_REDUCED_FIBER_{target_fiber}", combined_intensity_array)
            setattr(ad[0], f"OPTIMAL_REDUCED_ERR_{target_fiber}", combined_error_array)
            if symmetric_linefits:
                setattr(ad[0], f"WLS_SIMULTANEOUS_SYM_FIBER_{target_fiber}", combined_wavelength_array)
            else:
                setattr(ad[0], f"WLS_SIMULTANEOUS_FIBER_{target_fiber}", combined_wavelength_array)

            # Mark history
            fiber_list = ','.join(map(str, combine_fibers))
            gt.mark_history(ad, primname=self.myself(),
                        keyword='FIBER_COMBINATION',
                        comment=f"combined_fibers_{fiber_list}_to_{target_fiber}")
            
            log.info(f"Combined fibers {fiber_list} into fiber {target_fiber} for {ad.filename}")
            ad.update_filename(suffix=params["suffix"], strip=True)
        
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def barycentricCorrection(self, adinputs=None, **params):
        """
        Calculate barycentric velocity corrections for science observations.

        This primitive computes barycentric radial velocity corrections (BERV)
        for high-precision radial velocity measurements. It uses the barycorrpy
        library to calculate Earth's velocity projection toward the target at
        various exposure timestamps, accounting for exposure meter flux-weighted
        timing for maximum accuracy.

        The primitive calculates BERV at the exposure midpoint and flux-weighted
        midpoint using exposure meter data from both PC and FRD channels. It
        queries SIMBAD for target coordinates or uses telescope pointing directly.
        Timing corrections account for instrument-specific offsets between UTC
        filename timestamps and actual exposure times.

        Parameters
        ----------
        adinputs : list of AstroData
            Input AstroData objects containing 1D extracted science spectra with
            EXPOSUREMETER extension, timing information (UT_DATETIME, MJD), and
            telescope pointing data in headers.
        target_name : str, optional
            Target name substring for file filtering. Only files with OBJECT
            header matching this string are processed. If None, all files are
            processed. Default is None.
        simbad_target_name : str, optional
            SIMBAD-resolvable target name to override OBJECT header value. Use
            when the header target name differs from the SIMBAD catalog name.
            Only applies if target_name matches. Default is None.
        use_coords : bool, optional
            If True, use telescope pointing coordinates (TELRA, TELDEC) directly
            instead of querying SIMBAD for target coordinates. Recommended when
            target is not in SIMBAD or has unreliable proper motion data.
            Default is False.
        zp_pc : float, optional
            Zeropoint for PC (Precision Coupler) exposure meter channel in counts.
            If 0.0, automatically determined from median of 20 lowest readings in
            a 10-minute window around exposure. Units: counts. Default is 0.0.
        zp_frd : float, optional
            Zeropoint for FRD (Fiber Refractive Diffraction) exposure meter channel
            in counts. If 0.0, automatically determined from median of 20 lowest
            readings in a 10-minute window around exposure. Units: counts.
            Default is 0.0.
        start_time : str, optional
            Method for determining exposure start time. Options:
            - 'filename': Use UTC from filename (default)
            - 'mjd_start': Use telescope MJD at exposure start
            - 'mjd_end': Use telescope MJD at readout end minus exposure time
            Different methods account for varying instrument timing behaviors
            between red and blue arms. Default is 'filename'.
        suffix : str, optional
            Suffix to append to output filenames. Default is '_reduced'.

        Returns
        -------
        list of AstroData
            Modified AstroData objects with barycentric velocity corrections and
            exposure meter statistics added to the first extension header:

            BERV values (m/s):
            - BERV_MIDPOINT: BERV at nominal exposure midpoint
            - BERV_FLUXWEIGHTED_PC: BERV at flux-weighted midpoint (PC channel)
            - BERV_FLUXWEIGHTED_FRD: BERV at flux-weighted midpoint (FRD channel)
            - BERV_DIFFERENCE_PC: Difference between flux-weighted and nominal
            - BERV_DIFFERENCE_FRD: Difference between flux-weighted and nominal

            Timing information:
            - UTC_START: Corrected UTC start time (ISO format)
            - UTC_MIDPOINT: UTC midpoint time
            - UTC_FLUXWEIGHTED_PC: Flux-weighted UTC (PC channel)
            - UTC_FLUXWEIGHTED_FRD: Flux-weighted UTC (FRD channel)
            - UTC_CORRECTION: Applied time correction in seconds
            - JD_UTC_START: Julian date at start
            - JD_UTC_MIDPOINT: Julian date at midpoint
            - JD_UTC_FLUXWEIGHTED_PC: Flux-weighted JD (PC channel)
            - JD_UTC_FLUXWEIGHTED_FRD: Flux-weighted JD (FRD channel)

            Exposure meter statistics (counts):
            - COUNTS_PC_MIN/MAX/MEDIAN/STD: PC channel statistics
            - COUNTS_FRD_MIN/MAX/MEDIAN/STD: FRD channel statistics
            - COUNTS_PC_ZP: Applied PC zeropoint
            - COUNTS_FRD_ZP: Applied FRD zeropoint
            - SCALEFACTOR: Ratio of FRD to PC median counts

            Target information:
            - BERV_SIMBAD_TARGET: Target name used for BERV calculation

        References
        ----------
        .. [1] barycorrpy: https://github.com/shbhuk/barycorrpy
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Get parameters
        target_name = params.get('target_name')
        simbad_target_name = params.get('simbad_target_name')
        use_coords = params.get('use_coords')
        zp_pc = params.get('zp_pc')
        zp_frd = params.get('zp_frd')
        start_time = params.get('start_time')
        
        for ad in adinputs:
            # Read target name from header
            target = ad.object()
            log.debug(f'{target_name}, {target}')
            # Handle user-supplied target name logic
            if target_name is not None:
                if target_name in target:
                    log.fullinfo(f'Selected target name: {target}')
                    if simbad_target_name and len(simbad_target_name) > 0:
                        log.fullinfo(f'Replaced with user supplied name: {simbad_target_name}')
                        target = simbad_target_name
                else:
                    log.warning(f'Skip file {ad.filename}')
                    continue

            # Query SIMBAD for target coordinates
            result_table = Simbad.query_object(target)
            if result_table is None:
                log.warning(f'Target {target} not recognized by SIMBAD')
                if use_coords:
                    log.warning(f'Will use telescope pointing as coordinates.')
                else:
                    log.warning(f'Skip file {ad.filename} - consider using --use_coords option')
                    continue

            # Calculate time correction based on selected method
            utc_start = Time(ad.ut_datetime(), format='datetime', scale='utc')
            time_correction = TimeDelta(0, format='sec')

            mjd = ad[0].telescope_mjd(pretty=True)
            exptime = ad[0].exposure_time(pretty=True)
            if start_time == 'mjd_start':
                time_correction = mjd - utc_start
            elif start_time == 'mjd_end':
                if 'BLUE' in ad.tags:
                    # 52 and 100 come from legacy code, no docs about them
                    offset = exptime + TimeDelta(52.0, format='sec')
                else:
                    offset = exptime + TimeDelta(100.0, format='sec')
                time_correction = (mjd - offset) - utc_start

            # Extract times and exposure meter readings for a given exposure
            emeter = self._exposuremeterStats(ad, utc_start, exptime, zp_pc=zp_pc, zp_frd=zp_frd)
            log.fullinfo(f'Exposuremeter stats (PC channel): {emeter["pc"]["stats"]}')
            log.fullinfo(f'Exposuremeter stats (FRD channel): {emeter["frd"]["stats"]}')

            # Calculate BERV at start, mid, end of exposure
            utc_start += time_correction
            utc_mid = utc_start + exptime / 2
            utc_end = utc_start + exptime

            
            # Call barycorrpy for barycentric calculations
            barycorrpy_kwargs = {
                'lat': 19.823801,
                'longi': -155.469047,
                'alt': 4213,
                'ephemeris': 'de430',
                'zmeas': 0.0,
                'predictive': False,
                'leap_update': False
            }    
            if use_coords:
                ra  = ad[0].hdr.get('MAROONX TELESCOPE TELRA') * 15.0
                dec = ad[0].hdr.get('MAROONX TELESCOPE TELDEC')
                barycorrpy_kwargs['ra'] = ra
                barycorrpy_kwargs['dec'] = dec
                log.fullinfo(f'Using RA: {ra:.4f} deg, DEC: {dec:.4f} deg')
            else:
                barycorrpy_kwargs['starname'] = target
                log.fullinfo(f'Using target name: {target}')

            # BVC for nominal exposure midtime
            result1, warning, status = get_BC_vel(
                JDUTC=[utc_start.jd, utc_mid.jd, utc_end.jd],
                **barycorrpy_kwargs)

            # Average BVC over exposure weighted by exposure meter readings - BEST
            if emeter['pc']['times'] is not None:
                result3, JDUTCMID_pc, warning3, status3 = exposure_meter_BC_vel(
                    JDUTC=emeter['pc']['times'].jd + time_correction.jd,
                    expmeterflux=emeter['pc']['readings'],
                    **barycorrpy_kwargs)
            if emeter['frd']['times'] is not None:
                result4, JDUTCMID_frd, warning4, status4 = exposure_meter_BC_vel(
                    JDUTC=emeter['frd']['times'].jd + time_correction.jd,
                    expmeterflux=emeter['frd']['readings'],
                    **barycorrpy_kwargs)  

            # dvdt values in m/s/s
            m_s = units.Unit('m/s')
            m_s_s = units.Unit('m/s/s')
            BC_dvdt = [
                (result1[2] - result1[0]) * m_s / (exptime),
                (result1[2] - result1[1]) * m_s / (exptime/2),
                (result1[1] - result1[0]) * m_s / (exptime/2)]


            log.fullinfo(f'File: {ad.filename}')
            log.fullinfo(f'Target: {target}')
            log.fullinfo(f'Exptime: {exptime} sec')
            log.fullinfo(f'BERV dv/dt: '
                f'{BC_dvdt[0].to(m_s_s).value:.4f}, '
                f'{BC_dvdt[1].to(m_s_s).value:.4f}, ' 
                f'{BC_dvdt[2].to(m_s_s).value:.4f} m/s/s')

            # Mid values
            if emeter['pc']['times'] is not None:
                utc_fluxmid_pc = Time(JDUTCMID_pc, format="jd", scale="utc")
                utc_fluxmid_frd = Time(JDUTCMID_frd, format="jd", scale="utc")
                BC_fluxmid_pc = result3
                BC_fluxmid_frd = result4

            # If times are None, there was a gap in the photometer data and midpoints 
            # and BERVs are taken from the nominal midpoint.
            if emeter['pc']['times'] is None:
                utc_fluxmid_pc = utc_mid
                BC_fluxmid_pc = result1[1]

            if emeter['frd']['times'] is None:
                utc_fluxmid_frd = utc_mid
                BC_fluxmid_frd = result1[1]

            # Scale factor
            scale_factor = emeter['frd']['stats']['median'] / emeter['pc']['stats']['median']

            # Save header entries
            ad[0].hdr.set('UTC_START', f'{utc_start.isot}', 'xxxxx')
            ad[0].hdr.set('UTC_CORRECTION', f'{time_correction.sec:.1f}', 'xxxxx')
            ad[0].hdr.set('UTC_MIDPOINT', f'{utc_mid.isot}', 'xxxxx')
            ad[0].hdr.set('UTC_FLUXWEIGHTED_PC', f'{utc_fluxmid_pc.isot}', 'xxxxx')
            ad[0].hdr.set('UTC_FLUXWEIGHTED_FRD', f'{utc_fluxmid_frd.isot}', 'xxxxx')
            ad[0].hdr.set('JD_UTC_START', f'{utc_start.jd:.7f}', 'xxxxx')
            ad[0].hdr.set('JD_UTC_MIDPOINT', f'{utc_mid.jd:.7f}', 'xxxxx')
            ad[0].hdr.set('JD_UTC_FLUXWEIGHTED_PC', f'{utc_fluxmid_pc.jd:.7f}', 'xxxxx')
            ad[0].hdr.set('JD_UTC_FLUXWEIGHTED_FRD', f'{utc_fluxmid_frd.jd:.7f}', 'xxxxx')

            ad[0].hdr.set('BERV_SIMBAD_TARGET', target, 'xxxxx')
            ad[0].hdr.set('BERV_MIDPOINT', f'{result1[1]:.2f}', 'xxxxx')
            ad[0].hdr.set('BERV_FLUXWEIGHTED_PC', f'{BC_fluxmid_pc:.2f}', 'xxxxx')
            ad[0].hdr.set('BERV_FLUXWEIGHTED_FRD', f'{BC_fluxmid_frd:.2f}', 'xxxxx')
            ad[0].hdr.set('BERV_DIFFERENCE_PC', f'{(BC_fluxmid_pc - result1[1]):.2f}', 'xxxxx')
            ad[0].hdr.set('BERV_DIFFERENCE_FRD', f'{(BC_fluxmid_frd - result1[1]):.2f}', 'xxxxx')
            # ad[0].hdr.set('BERV_DVDT', f'{np.mean(BC_dvdt):.2f}', 'xxxxx')

            ad[0].hdr.set('COUNTS_PC_MIN', f'{emeter['pc']['stats']['min']:.2f}', 'xxxxx')
            ad[0].hdr.set('COUNTS_PC_MAX', f'{emeter['pc']['stats']['max']:.2f}', 'xxxxx')
            ad[0].hdr.set('COUNTS_PC_MEDIAN', f'{emeter['pc']['stats']['median']:.2f}', 'xxxxx')
            ad[0].hdr.set('COUNTS_PC_STD', f'{emeter['pc']['stats']['std']:.2f}', 'xxxxx')
            ad[0].hdr.set('COUNTS_PC_ZP', f'{emeter['pc']['stats']['zeropoint']:.2f}', 'xxxxx')

            ad[0].hdr.set('COUNTS_FRD_MIN', f'{emeter['frd']['stats']['min']:.2f}', 'xxxxx')
            ad[0].hdr.set('COUNTS_FRD_MAX', f'{emeter['frd']['stats']['max']:.2f}', 'xxxxx')
            ad[0].hdr.set('COUNTS_FRD_MEDIAN', f'{emeter['frd']['stats']['median']:.2f}', 'xxxxx')
            ad[0].hdr.set('COUNTS_FRD_STD', f'{emeter['frd']['stats']['std']:.2f}', 'xxxxx')
            ad[0].hdr.set('COUNTS_FRD_ZP', f'{emeter['frd']['stats']['zeropoint']:.2f}', 'xxxxx')
            ad[0].hdr.set('SCALEFACTOR', f'{scale_factor:.1f}', 'xxxxx')

            log.fullinfo(f"Barycentric velocity calculation completed for {ad.filename}")
            ad.update_filename(suffix=params["suffix"], strip=True)

        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        log.debug(gt.log_message("primitive", self.myself(), "complete"))
        return adinputs

    def _exposuremeterStats(self, ad, utc_start, exptime, zp_pc=0.0, zp_frd=0.0):
        """
        Extract and process exposure meter readings.
        
        This function retrieves exposure meter data from a pandas dataframe, applies
        zeropoint corrections, performs outlier filtering, and calculates statistics
        for both PC and FRD channels.
        
        Parameters
        ----------
        ad : AstroData object
            Input AstroData object containing EXPOSUREMETER extension.
        utc_start : astropy.time.Time
            UTC start time of the exposure.
        exptime : float or astropy.time.TimeDelta
            Exposure time in seconds.
        zp_pc : float, optional
            Zeropoint for PC channel. If 0.0 (default), automatically determined
            from the median of the 20 lowest values in a 10-minute window around
            the exposure.
        zp_frd : float, optional
            Zeropoint for FRD channel. If 0.0 (default), automatically determined
            from the median of the 20 lowest values in a 10-minute window around
            the exposure.
        
        Returns
        -------
        dict
            Nested dictionary containing processed exposure meter data:
            
            - 'pc'/'frd' : dict
                Channel-specific data containing:
                - 'times' : astropy.time.Time or None
                    Timestamps for readings during exposure
                - 'readings' : numpy.ndarray or nan
                    Zeropoint-corrected flux readings (negative values set to 0)
                - 'stats' : dict
                    Statistical measures:
                    - 'min', 'max', 'median', 'std' : float
                        Basic statistics of the readings
                    - 'zeropoint' : float
                        Applied zeropoint value
        """
        log = self.log
        exposuremeter = ad.EXPOSUREMETER.to_pandas(index='Timestamp')  
        exposuremeter.index = pd.to_datetime(exposuremeter.index)
        
        # Reference zeropints
        ref_zp_pc = ad.EXPOSUREMETER.meta['header']['ZP_PC']
        ref_zp_frd = ad.EXPOSUREMETER.meta['header']['ZP_FRD']

        # 5 min window around exposure for auto-zp determination
        dt1 = TimeDelta(exptime, format='sec')
        dt3 = TimeDelta(300, format='sec')
        utc_end = utc_start + dt1
        
        # Auto-determine zeropoints if not provided
        n_cutoff = 20
        start_cut = (utc_start - dt3).iso
        end_cut = (utc_end + dt3).iso
        if zp_frd == 0.0:
            flux_frd = exposuremeter.loc[start_cut:end_cut]['Flux FRD Channel']
            zp_frd = np.nanmedian(np.sort(flux_frd)[:n_cutoff])
            log.fullinfo(f'Automatic zeropoint determination for FRD channel: {zp_frd}')
        if zp_pc == 0.0:
            flux_pc = exposuremeter.loc[start_cut:end_cut]['Flux PC Channel']
            zp_pc = np.nanmedian(np.sort(flux_pc)[:n_cutoff])
            log.fullinfo(f'Automatic zeropoint determination for PC channel: {zp_pc}')
        
        # Check if zp are valid
        if np.isnan(zp_frd):
            log.warning("Automatic zeropoint determination for FRD channel failed.")
            zp_frd = 0.0
        if np.isnan(zp_pc):
            log.warning("Automatic zeropoint determination for PC channel failed.")
            zp_pc = 0.0
        if abs(zp_frd - ref_zp_frd) > 0.2 * ref_zp_frd:
            log.warning(f"Automatic zeropoint determination for FRD, {zp_frd}, "
                        f"is 20% off of reference value of {ref_zp_frd}.")
        if abs(zp_pc - ref_zp_pc) > 0.2 * ref_zp_pc:
            log.warning(f"Automatic zeropoint determination for FRD, {zp_pc}, "
                        f"is 20% off of reference value of {ref_zp_pc}.")

        # Extract readings during exposure
        result_pc = exposuremeter.loc[utc_start.iso:utc_end.iso]['Flux PC Channel']
        result_frd = exposuremeter.loc[utc_start.iso:utc_end.iso]['Flux FRD Channel']
        number_pc = result_pc.count()
        number_frd = result_frd.count()

        # Apply zp correction and outlier filtering
        times_pc = Time(result_pc.index.values, format='datetime64', scale='utc')
        readings_pc = result_pc.values.flatten() - zp_pc
        median_pc = medfilt(readings_pc, 3)
        outlier = np.where(np.abs(readings_pc - median_pc) / median_pc > 2)
        readings_pc[outlier] = median_pc[outlier]        
        if np.sum(outlier) > 0:
            log.warning(f'Replaced {len(outlier)} outlier value(s) in PC dataset')        
        readings_pc[readings_pc < 0] = 0.0
        
        times_frd = Time(result_frd.index.values, format='datetime64', scale='utc')
        readings_frd = result_frd.values.flatten() - zp_frd
        median_frd = medfilt(readings_frd, 3)
        outlier = np.where(np.abs(readings_frd - median_pc) / median_frd > 2)
        readings_frd[outlier] = median_frd[outlier]
        if np.sum(outlier) > 0:
            log.warning(f'Replaced {len(outlier)} outlier value(s) in FRD dataset')
        readings_frd[readings_frd < 0] = 0.0

        # Check for large time gaps in data
        if number_frd > 2 and number_pc > 2:
            gaps = [
                (times_frd[0] - utc_start).sec, 
                (utc_end - times_frd[-1]).sec,
                result_frd.index.to_series().diff().dt.total_seconds().fillna(0).max()
                ]
            maxgap = max(gaps)
            if maxgap > 30:
                log.warning(f'{maxgap:.1f} sec gap found in exposuremeter data. Photometric calculations abandoned.')
            elif maxgap > 10:
                logger.warning(f'{maxgap:.1f} sec gap found in exposuremeter data.')
        else:
            times_pc = None
            times_frd = None
            readings_pc = np.nan
            readings_frd = np.nan

        # Calculate statistics
        emeter_stats = {
            'pc': {
                'times': times_pc,
                'readings': readings_pc,
                'stats': {
                    'min': np.nanmin(readings_pc),
                    'max': np.nanmax(readings_pc),
                    'median': np.nanmedian(readings_pc),
                    'std': np.nanstd(readings_pc),
                    'zeropoint': zp_pc
                }
            },
            'frd': {
                'times': times_frd,
                'readings': readings_frd,
                'stats': {
                    'min': np.nanmin(readings_frd),
                    'max': np.nanmax(readings_frd),
                    'median': np.nanmedian(readings_frd),
                    'std': np.nanstd(readings_frd),
                    'zeropoint': zp_frd
                }
            }
        }
        return emeter_stats


##############################################################################
# Below are the helper functions for the primitives in this module           #
##############################################################################

def _make_b_spline_from_pars(p, kind=5):
    disp_params = []
    for par in p:
        if "disp_" in par:
            disp_params.append(p[par].value)
    disp_params = np.array(disp_params, dtype=float)
    t = []
    for par in p:
        if "knot_" in par:
            t.append(p[par].value)
    knots = np.array(t, dtype=float)

    return scipy.interpolate.BSpline(knots, disp_params, kind)

def _peak_to_wavelength(m, pars):
    return (2.0 * (pars["l"]) * np.cos(pars["theta"]) * pars["n"] / m) * 1e6

def _peak_to_wavelength_spline(mm, pars):
    spl = _make_b_spline_from_pars(pars)
    return (
        2.0 * (pars["l"]  - spl(1 / mm)*pars["l"]) * np.cos(pars["theta"]) * pars["n"] / mm
    ) * 1e6

def _fc2min(p, m, etalonwl):
    # residuals are in 'nm' not m/s. Good? bad? Should we normalize?
    return _peak_to_wavelength_spline(m, p) - etalonwl