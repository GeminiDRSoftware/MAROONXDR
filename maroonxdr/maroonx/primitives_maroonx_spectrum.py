"""
This module contains primitives to generate wavelength
calibration solutions from reduced 1-D spectra.
"""
# ------------------------------------------------------------------------------
import multiprocessing
import time
import traceback

from astropy.table import Table
import numpy as np
import pandas as pd
import scipy

from geminidr.core import Spect
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

from . import parameters_maroonx_spectrum
from .maroonx_fit import maroonx_fit
from . import maroonx_utils
from .primitives_maroonx_echelle import MAROONXEchelle
from .maroonx_echellespectrum.maroonxspectrum import MXSpectrum
from .maroonx_echellespectrum.wavelengthsolution import WavelengthSolution
# ------------------------------------------------------------------------------

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
            # See which fibers and orders we are extracting
            if fibers is None:
                log.warning("No fibers specified.  Finding all etalon fibers and extracting those.")
                fibers_list = []
                # Read the HIERARCH FIBER keywords from the phu
                if ad.phu['FIBER1'] == 'Etalon':
                    fibers_list.append(1)
                if ad.phu['FIBER2'] == 'Etalon':
                    fibers_list.append(2)
                if ad.phu['FIBER3'] == 'Etalon':
                    fibers_list.append(3)
                if ad.phu['FIBER4'] == 'Etalon':
                    fibers_list.append(4)
                if ad.phu['FIBER5'] == 'Etalon':
                    fibers_list.append(5)
                fibers = tuple(fibers_list)

            if orders is None:
                log.warning("No orders specified.  Extracting all orders.")

            # If multithreading is on, disable plotting
            if multithreading and show_plots:
                log.warning("Multithreading is on.  Disabling plotting.")
                show_plots = False

            errors = []
            results = []
            log.fullinfo(f"Extracting fibers {fibers}")

            # Create worker pools for multithreading
            if multithreading:
                pool = multiprocessing.get_context("fork").Pool()
                log.fullinfo(f"Using {pool._processes} processes for extraction")
                for fiber, order, data, guess in maroonx_utils.load_recordings(
                ad, guess_file, fibers, orders):
                    log.fullinfo(f"Fitting fiber {fiber} order {order}")

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
                        log.warning(f"Error extracting fiber {fiber} order {order}: {exc}")

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
                    log.fullinfo(f"Fitting fiber {fiber} order {order}")

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

                    # Run iterative fit using the data yielded by the generator function in serial
                    output = maroonx_fit.iterative_fit(
                            input_spectrum = data,
                            degree_sigma = degree_sigma,
                            degree_width = degree_width,
                            iterations = iterations,
                            guess_spectrum = guess,
                            fiber="{}_{}".format(fiber, order),
                            plot_path=plot_path,
                            use_sigma_lr = use_sigma_lr,
                            show_plots = show_plots)

                    # Save the results in a manner anologoous to the save_results callback function
                    output.fiber = fiber
                    output.order = order
                    output.recording_time = 0.0  # TODO extract from data set
                    results.append(output)

            # Record the time taken
            end_time = time.time()
            log.fullinfo(f"Finished extracting etalon lines in {end_time - start_time:.2f} seconds")

            if results is not []:
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
        Computes a new spline-based dynamical wavelength solution for each fiber 
        using etalon parameters from a config file and an initial wavelength solution.
        
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
            # Load the etalon spectrum
            mx_spectrum = MXSpectrum(ad, etalon_peaks_symmetric=symmetric_linefits)
            log.info(f'Processing etalon file: {ad.filename}')
            
            # Determine fibers to process if not provided
            if fibers is None:
                fibers = []
                for fiber_num, fiber_type in enumerate(ad.fiber_setup(), start=1):
                    if fiber_type == 'Etalon':
                        fibers.append(fiber_num)
                log.info(f'Found etalon fibers: {fibers}')
            
            if ref_file is not None:
                # If a reference file and reference fiber is given, offset the etalon
                # position on all fibers by the offset found in the reference fiber.
                # ref_fiber = params.get('ref_fiber')
                # TODO: Ask Andreas what this is supposed to do because currently
                # we do not know how this works in the old pipeline
                raise NotImplementedError("Reference file not implemented yet")

            # Load reference wavelength solution from the config
            refwavelength_file = maroonx_utils.get_refwavelength_filename(ad)
            log.info(f'Loading reference wavelength file: {refwavelength_file}')
            
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
                    log.info(f'Apply static wavelength solution for peak number identification in fiber {fiber}.')
                    mx_spectrum.spectra[fiber].apply_wavelength_vector()

            log.info(f'Apply etalon parameters from file {refwavelength_file}.')
            parameters = maroonx_utils.load_params_from_fits(refwavelength_file, ext_name='PARAMETERS')
            for fiber in fibers:
                mx_spectrum.spectra[fiber].etalon_pars = parameters

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

                if fiber == 5:
                    residuals = _fc2min(parameters, peak_data["M"].values, peak_data["WAVELENGTH_BY_THAR"].values) / peak_data[
                        "WAVELENGTH_BY_THAR"].values * 3e8
                    bad = np.where(np.abs(residuals-np.nanmedian(residuals)) > 4.0 * np.nanstd(residuals))
                    residuals[bad] = np.nan
                    inst_drift = np.nanmean(residuals)
                    drifts[fiber] = inst_drift
                else:
                    residuals = _fc2min(parameters, peak_data["M"].values, peak_data["WAVELENGTH_BY_THAR"].values) / peak_data[
                        "WAVELENGTH_BY_THAR"].values * 3e8
                    bad = np.where(np.abs(residuals-np.nanmedian(residuals)) > 4.0 * np.nanstd(residuals))
                    residuals[bad] = np.nan
                    drift = np.nanmean(residuals)
                    drifts[fiber] = drift
                    #drifts = np.append(drifts, drift)
                
            # ============================================================== OK line


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
                for o in spectra.orders:

                    order_mask = spectra.peak_data["ORDER"] == o
                    center = spectra.peak_data["CENTER"][order_mask]

                    if center.values[0] < center.values[10]:
                        x = (center.values)
                        y = (_peak_to_wavelength_spline(spectra.peak_data["M"][order_mask],
                                                    spectra.etalon_pars).values)
                    else:
                        x = (center.values)[::-1]
                        y = (_peak_to_wavelength_spline(spectra.peak_data["M"][order_mask],
                                                    spectra.etalon_pars).values)[::-1]
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

                new_peak_data.append(spectra.peak_data)

                # Concatenate wavelengths arrays for each fiber and save them
                wave_arrays = [wave[fiber][o] for o in wave[fiber].keys()]
                setattr(ad[0], f'WLS_DYNAMIC_FIBER_{fiber}', np.stack(wave_arrays))

            for fiber in {1, 2, 3, 4, 5} - set(fibers):
                # If fiber is not in the list of fibers, save an empty array
                setattr(ad[0], f'WLS_DYNAMIC_FIBER_{fiber}', np.zeros((1, 1)))

            # Collect and reformat updated peak data
            new_peak_data = pd.concat(new_peak_data, ignore_index=True)
            new_peak_data.sort_values(by=['FIBER', 'ORDER', 'M'], inplace=True)
            ad[0].NEW_PEAKS = Table.from_pandas(new_peak_data)

            # Save meassured drift in header entries
            for fiber in fibers:
                ad[0].hdr[f'DRIFT_FIBER_{fiber}'] = (drifts[fiber], "Drift in m/s")
                log.info(f"Drift for fiber {fiber}: {drifts[fiber]} m/s")

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

# ============================================================================================================================
# ============================================================================================================================

def fitAndApplyEtalonWlsOLD(self, adinputs, plot_path=None, ref_file=None, ref_fiber=None, symmetric_linefits=False):
    """
    This step computes a new spline-based dynamical wavelength solution for each fiber using etalon paramters that
    are provided as a config file, and a initial solution generated from a DTTTE file.  The lines in the spectra
    are identified based on the initial spectrum, which assigns a corresponding wavelength to each pixel.  Then
    this is fitted to a 30 knot spline to find the new dynamic wavelength solution.

    Parameters:
    -----------
    adinputs: list of AstroData objects with 1D box extracted spectra (PEAKS and POLY extensions)
    fibers: list of ints
        Fibers containing Etalon spectra
    plot_path: str
        If path is not none, save plots to this path.
    ref_file: str
        Absolute path and filename containing reduced and fitted etalon spectra.
        This input is only used for the drift correction step.
    ref_fiber: int
        Which fiber to use as the reference fiber if a reference spectrum was provided.
        This input is only used for the drift correction step.

    Returns:
    --------
    adinputs with a new extension containing the spline based wavelength solution
    """
    log = self.log
    log.debug(gt.log_message("primitive", self.myself(), "starting"))

    # Determine tag for log and plot filename
    if symmetric_linefits:
        tag  = '_spline_symmetrical'
    else:
        tag = '_spline'
    if ref_file is not None:
        tag = tag + '_ref'

    for adinput in adinputs:
        # Create pdf for plots
        filename = adinput.filename
        pp = PdfPages(filename[:-4] + tag + '.pdf')
        try:
            # Load etalon spectrum
            mx_obj = MXSpectrum(adinput, etalon_peaks_symmetric=symmetric_linefits)

            log.info(f'Etalon file: {filename}')
            fiber2_obj = mx_obj.spectra[2]
            fiber3_obj = mx_obj.spectra[3]
            fiber4_obj = mx_obj.spectra[4]
            fiber5_obj = mx_obj.spectra[5]
            
            fiber2_peak_data = fiber2_obj.peak_data
            print(fiber2_peak_data)
        except Exception as e:
            log.error(f'Error processing file: {filename}')
            log.error(f'Exception: {e}')
            continue

    if ref_file is not None:
        ref_peaks = ref_file[0].PEAKS[ref_file[0].PEAKS['FIBER'] == ref_fiber]
        peak_centers = ref_peaks["CENTER"]
        # TODO: Ask Andreas what this is supposed to do because currently
        # we do not know how this works in the old pipeline

    for ad in adinputs:
        # Set the location of the plot file
        if plot_path is not None:
            plot_file = plot_path + '/' + adinputs[0].filename + tag + '.png'
            
        # Load the reference spectrum from config
        wavelength_file = maroonx_utils.get_refwavelength_filename(ad)
        wave_dict = ad.open(wavelength_file)
        log.fullinfo(f"Using reference file {wavelength_file} for wavelength solution")
        if symmetric_linefits:
            etalon = MXSpectrum(ad, etalon_peaks_symmetric=True)
        else:
            etalon = MXSpectrum(ad)
        for fiber in fibers:
            # Guess the peak numbers in the measured spectrum
            peak_numbers = guess_peak_numbers(self, reduced_fiber, peak_data, poly_data)
        return adinputs

