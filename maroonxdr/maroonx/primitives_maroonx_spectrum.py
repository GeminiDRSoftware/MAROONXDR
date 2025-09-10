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
import numpy as np
import pandas as pd
import scipy
from scipy.interpolate import LSQUnivariateSpline, interp1d

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
        Loads the static wavelength solution from the config file.

        The config file is located in the lookup directory, see wavelengthdb.py
        
        Parameters:
        -----------
        adinputs: list 
            AstroData objects with 1D box extracted spectra
        fibers: list or tuple
            Fibers containing Etalon spectra to process

        Returns:
        --------
        adinputs with extensions containing the static wavelengths solution
        """

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # requested fibers
        fibers = params['fibers']
        
        all_fibers = [1, 2, 3, 4, 5]
        if fibers is None:
            fibers = all_fibers

        for ad in adinputs:
            # Load static wavelength solution from the config
            statwavelength_file = maroonx_utils.get_statwavelength_filename(ad)
            log.info(f'Loading static wavelength file: {statwavelength_file}')

            for fiber in all_fibers:
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

        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def getPeaksAndPolynomials(self, adinputs=None, guess_file=None,
                            fibers=(), orders=(), degree_sigma=4, degree_width=2,
                            use_sigma_lr=True, show_plots=False,
                            plot_path="", multithreading=False,
                            iterations = 5, **params):
        """
        Extracts the etalon positions from the 1D spectra, determines the centroid,
        and fits polynomials.  This is done by finding the peak and fitting them using a
        box convolved with 2 Gaussians  The sigmas of the Gaussians are fitted to a low
        order polynomial, along with the width of the box.

        Below information taken from the MAROON-X Data Handbook:
        https://sites.google.com/uchicago.edu/maroonx-data-handbook/data-reduction/wavelength-calibration/etalon-line-fitting

        The line profiles are modeled using a box
        convolved with 2 Gaussians on either side (sum of 2 erf functions).
        In order to remove degeneracies between the width of the box (representing
        the width of the entance slit) and the sigma/FWHM of the Gaussians (representing
        the field and wavelength dependent abberations), a composite model for each fiber/order
        is built from all the individual etalon lines in a fiber/order plus a smooth background.
        Only the line intensities and positions are fitted individually for each  etalon line.
        Box widths, Gaussian sigmas, and the background are constrained to a low order polynomial
        across a given fiber/ order.  This is motivated by the fact that these quantities only vary
        slowly and steadily along an order.

        Parameters:
        -----------
        adinputs: list of AstroData objects with 1D box extracted spectra
        guess_file: str
            Name of file containing initial guess spectrum.
        fibers: list of ints
            Fibers to fit
        orders: list of ints
            Orders to fit
        degree_sigma: int
            Degree of the sigma polynomial
        degree_width: int
            Degree of the width polynomial
        use_sigma_lr: bool
            Use different polynomial degrees for left and right side of wings
        show_plots: bool
            Show plots of the etalon line fits to be used for debugging.
        plot_path: str
            If show_plots is True, save plots to this path.
        multithreading: bool
            Use multithreading to speed up extraction.  Disables plotting.
        iterations: int
            Maximum number of iterations on the fit.

        Returns:
        --------
        adinputs with extensions containing the peak parameters and polynomial parameters
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

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
        Calculates spline-based dynamic wavelength solutions from etalon spectra.
        This primitive processes etalon files to create reference wavelength solutions
        needed for drift correction in science frames.
        
        Parameters:
        -----------
        adinputs: list
            AstroData objects with 1D extracted spectra (PEAKS and POLY extensions)
        fibers: list or tuple
            Fibers containing Etalon spectra to process
        symmetric_linefits: bool
            Whether to use symmetric line fitting
        n_knots: int
            Number of knots for the spline fit (default: 30)
        thar: bool
            Whether to apply ThAr wavelength solution to etalon frame
        ref_file: str
            Reference file with previously reduced etalon spectra (for drift correction)
        ref_fiber: int
            Fiber to use as reference when ref_file is provided
            
        Returns:
        --------
        adinputs: list
            Same AstroData objects with new wavelength solution
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

        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def applyWavelengthSolution(self, adinputs=None, **params):
        """
        Applies drift-corrected wavelength solutions to science spectra using 
        simultaneous etalon measurements. This primitive uses pre-calculated 
        etalon wavelength solutions to correct for instrumental drifts.
        
        Parameters:
        -----------
        adinputs: list
            AstroData objects with 1D extracted science spectra (PEAKS and POLY extensions)
        fibers: list or tuple
            Science fibers to process (default: [2,3,4])
        ref_fiber: int
            Fiber containing simultaneous etalon spectrum for drift measurement (default: 5)
        symmetric_linefits: bool
            Whether to use symmetric line fitting (default: False)
        n_knots: int
            Number of knots for the spline fit (default: 30)
        etalon_file: str
            Path to corresponding etalon file with dynamic wavelength solutions
            
        Returns:
        --------
        adinputs: list
            Same AstroData objects with drift-corrected wavelength solutions 
            stored as WLS_SCIENCE_FIBER_X extensions
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
        Combines data from multiple science fibers into one spectrum using weighted averaging.
        
        This primitive scales and interpolates fiber data onto a common wavelength grid,
        performs kappa-sigma clipping to reject outliers, and combines using inverse
        variance weighting.
        
        Parameters
        ----------
        adinputs: AstroData object(s) with extracted fiber data
        combine_fibers: (list) list of fiber numbers to combine
        symmetric_linefits: (bool) use symmetrical line fits (affects wavelength source)
        kappa_sigma: (float) sigma clipping threshold for outlier rejection
        max_clips: (int) maximum pixels to clip per order before increasing kappa_sigma
        
        Returns
        -------
        adinputs with combined fiber data in new extension with {target_fiber} suffix.
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
                setattr(ad[0], f"WAVELENGTH_SYM_FIBER_{target_fiber}", combined_wavelength_array)
            else:
                setattr(ad[0], f"WAVELENGTH_FIBER_{target_fiber}", combined_wavelength_array)

            # Mark history
            fiber_list = ','.join(map(str, combine_fibers))
            gt.mark_history(ad, primname=self.myself(),
                        keyword='FIBER_COMBINATION',
                        comment=f"combined_fibers_{fiber_list}_to_{target_fiber}")
            
            log.info(f"Combined fibers {fiber_list} into fiber {target_fiber} for {ad.filename}")
            ad.update_filename(suffix=params["suffix"], strip=False)
        
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

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

