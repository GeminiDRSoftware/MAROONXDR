"""
This module contains primitives to generate wavelength
calibration solutions from reduced 1-D spectra.
"""
# ------------------------------------------------------------------------------
import multiprocessing
import time
import traceback

from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from geminidr.core import Spect
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

from . import parameters_maroonx_spectrum
from .maroonx_fit import maroonx_fit
from . import maroonx_utils
from .primitives_maroonx_echelle import MAROONXEchelle
from .maroonx_echellespectrum.maroonxspectrum import MXSpectrum
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
        adinputs: list of AstroData objects with 1D box extracted spectra

        Returns:
        --------
        adinputs with extensions containing the static wavelengths solution
        """

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        fibers = params['fibers']
        if fibers is None:
            fibers = [1, 2, 3, 4, 5]

        for ad in adinputs:
            # Load static wavelength solution from the config
            statwavelength_file = maroonx_utils.get_statwavelength_filename(ad)
            log.info(f'Loading static wavelength file: {statwavelength_file}')

            for fiber in fibers:
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
        return

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

    def fitAndApplyEtalonWls(self, adinputs, plot_path=None, ref_file=None, ref_fiber=None, symmetric_linefits=False):
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


# ============================================================================================================================
# ============================================================================================================================

def fitAndApplyEtalonWlsNEW(self, adinputs=None, fibers=None, 
                         ref_file=None, ref_fiber=None, thar=False,
                         symmetric_linefits=False, n_knots=30, n_sigma_clip=4.0):
    """
    Computes a new spline-based dynamical wavelength solution for each fiber 
    using etalon parameters from a config file and an initial wavelength solution.
    
    Parameters:
    -----------
    adinputs: list
        AstroData objects with 1D extracted spectra (PEAKS and POLY extensions)
    fibers: list or tuple
        Fibers containing Etalon spectra to process
    ref_file: str
        Reference file with previously reduced etalon spectra (for drift correction)
    ref_fiber: int
        Fiber to use as reference when ref_file is provided
    thar: bool
        Whether to apply ThAr wavelength solution to etalon frame
    symmetric_linefits: bool
        Whether to use symmetric line fitting
    n_knots: int
        Number of knots for the spline fit (default: 30)
    n_sigma_clip: float
        Sigma clipping threshold for outlier rejection (default: 4.0)
        
    Returns:
    --------
    adinputs: list
        Same AstroData objects with new wavelength solution
    """
    log = self.log
    log.debug(gt.log_message("primitive", self.myself(), "starting"))
    
    # Process each input file
    for ad in adinputs:
        
        # Load the etalon spectrum
        mx_spectrum = MXSpectrum(ad, etalon_peaks_symmetric=symmetric_linefits)
        log.info(f'Processing etalon file: {ad.filename}')
        
        # Determine fibers to process if not provided
        if fibers is None:
            fibers = []
            for fiber_num, fiber_type in enumerate(ad.fiber_setup()):
                if fiber_type == 'Etalon':
                    fibers.append(fiber_num)
            log.info(f'Found etalon fibers: {fibers}')
        
        if ref_file is not None:
            # ref_peaks = ref_file[0].PEAKS[ref_file[0].PEAKS['FIBER'] == ref_fiber]
            # peak_centers = ref_peaks["CENTER"]
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

        # Create wavelength solutions dict to store results
        wavelength_solutions = {}

        # Process each fiber
        for fiber in fibers:
            etalon_spec = mx_spectrum.spectra[fiber]
            log.info(f'Processing fiber {fiber}')
            
            # Load the reference wavelength solution for this fiber
            ref_wls = None
            try:
                # Access the appropriate extension from the reference file
                from astropy.io import fits
                with fits.open(wavelength_file) as hdul:
                    ext_name = f'FIBER_{fiber}'
                    if ext_name in [hdu.name for hdu in hdul]:
                        ext = hdul[ext_name]
                        
                        # Create WavelengthSolution object
                        ref_wls = WavelengthSolution(
                            x_norm=ext.data['X_NORM'][0],
                            orders=ext.data['ORDERS'][0],
                            weights=ext.data['WEIGHTS'][0],
                            wavelengths=ext.data['WAVELEN'][0],
                            poly_deg_x=ext.header['POLYDEGX'],
                            poly_deg_y=ext.header['POLYDEGY'],
                            max_x=ext.header['MAXX']
                        )
                        log.info(f'Loaded reference wavelength solution for fiber {fiber}')
                    else:
                        log.warning(f'No wavelength solution found for fiber {fiber} in {wavelength_file}')
            except Exception as e:
                log.warning(f'Error loading reference wavelength solution for fiber {fiber}: {str(e)}')
            
            # Apply wavelength solution to spectrum
            if ref_wls is not None:
                log.info(f'Applying reference wavelength solution to fiber {fiber}')
                etalon_spec.apply_wavelength_solution(ref_wls)
            else:
                # Fall back to wavelength vector method
                log.info(f'Using wavelength vector for fiber {fiber}')
                etalon_spec.apply_wavelength_vector()
            
            # Guess etalon peak numbers
            log.info(f'Guessing etalon peak numbers for fiber {fiber}')
            peak_data = etalon_spec.guess_peak_numbers()
            
            # Create spline-based wavelength solutions for each order
            fiber_wavelengths = {}
            
            for order in etalon_spec.orders:
                try:
                    # Get peak data for this order
                    if order not in peak_data.index.levels[0]:
                        log.warning(f'No peak data for fiber {fiber} order {order}')
                        continue
                        
                    # Extract centers and wavelengths
                    centers = peak_data.loc[order, 'center'].values
                    wavelengths = peak_data.loc[order, 'wavelength'].values
                    
                    # Ensure data is sorted by center (important for spline fitting)
                    sort_idx = np.argsort(centers)
                    centers = centers[sort_idx]
                    wavelengths = wavelengths[sort_idx]
                    
                    # Create knots for spline
                    knots = np.linspace(np.min(centers) + 1, np.max(centers) - 1, n_knots)[1:-1]
                    
                    # Fit initial spline
                    from scipy.interpolate import LSQUnivariateSpline
                    spline = LSQUnivariateSpline(centers, wavelengths, knots, k=3)
                    
                    # Calculate residuals
                    residuals = wavelengths - spline(centers)
                    residuals_ms = residuals / wavelengths * 3.0e8  # Convert to m/s
                    
                    # Sigma clip outliers
                    std_residuals = np.nanstd(residuals_ms)
                    good_idx = np.abs(residuals_ms) < n_sigma_clip * std_residuals
                    
                    # Refit with outliers removed
                    if np.sum(good_idx) < len(centers):
                        log.info(f'Removed {np.sum(~good_idx)} outliers from fiber {fiber} order {order}')
                        spline = LSQUnivariateSpline(centers[good_idx], wavelengths[good_idx], knots, k=3)
                    
                    # Create full wavelength vector for all pixels
                    pixels = np.arange(len(etalon_spec.data.loc[order, 'wavelength']))
                    order_wavelengths = spline(pixels)
                    
                    # Handle extrapolation at edges if needed
                    # If spline returns values outside the expected range, use linear extrapolation
                    from scipy.interpolate import interp1d
                    
                    # Check if the first or last few values need extrapolation
                    edge_extrap_needed = False
                    if np.min(pixels) < np.min(centers) or np.max(pixels) > np.max(centers):
                        edge_extrap_needed = True
                        
                    if edge_extrap_needed:
                        # Create extrapolation function using the first/last few points
                        extrap_fn = interp1d(
                            np.concatenate([centers[:2], centers[-2:]]), 
                            np.concatenate([wavelengths[:2], wavelengths[-2:]]),
                            fill_value="extrapolate"
                        )
                        
                        # Apply extrapolation to edge pixels
                        left_mask = pixels < np.min(centers)
                        right_mask = pixels > np.max(centers)
                        if np.any(left_mask):
                            order_wavelengths[left_mask] = extrap_fn(pixels[left_mask])
                        if np.any(right_mask):
                            order_wavelengths[right_mask] = extrap_fn(pixels[right_mask])
                    
                    # Store wavelength solution for this order
                    fiber_wavelengths[str(order)] = order_wavelengths
                    log.info(f'Created wavelength solution for fiber {fiber} order {order}')
                    
                except Exception as e:
                    log.warning(f'Error creating wavelength solution for fiber {fiber} order {order}: {str(e)}')
            
            # Store all wavelengths for this fiber
            wavelength_solutions[f'fiber_{fiber}'] = fiber_wavelengths
        
        # Calculate any instrumental drifts
        drift_values = {}
        for fiber in fibers:
            try:
                etalon_spec = mx_spectrum.spectra[fiber]
                residuals = []
                
                # Collect all residuals for this fiber
                for order in etalon_spec.orders:
                    if order in peak_data.index.levels[0]:
                        m_values = peak_data.loc[order, 'm'].values
                        thar_wl = peak_data.loc[order, 'wavelength_by_thar'].values
                        
                        # Calculate residuals from model
                        res = (etalon_spec.peak_to_wavelength(m_values) - thar_wl) / thar_wl * 3.0e8
                        residuals.extend(res)
                
                if residuals:
                    # Filter out outliers before calculating mean drift
                    residuals = np.array(residuals)
                    good_idx = np.abs(residuals - np.nanmedian(residuals)) < 4.0 * np.nanstd(residuals)
                    drift = np.nanmean(residuals[good_idx])
                    drift_values[fiber] = drift
                    log.info(f'Calculated instrument drift for fiber {fiber}: {drift:.1f} m/s')
            except Exception as e:
                log.warning(f'Error calculating drift for fiber {fiber}: {str(e)}')
        
        # Save wavelength solutions and drift values to the AstroData object
        # Create wavelength extension
        from astropy.table import Table
        
        # Save wavelength solutions for each fiber
        for fiber, fiber_wavelengths in wavelength_solutions.items():
            fiber_num = int(fiber.split('_')[1])
            
            # Convert to table format
            wl_data = []
            for order, wavelengths in fiber_wavelengths.items():
                wl_data.append({'ORDER': int(order), 'WAVELENGTH': wavelengths})
            
            # Create table
            wl_table = Table(wl_data)
            
            # Add to AstroData
            extname = f'WAVELENGTH_FIBER_{fiber_num}'
            ad[0].append(wl_table, extname=extname)
            log.info(f'Added wavelength solution for fiber {fiber_num} to {ad.filename}')
        
        # Save drift values
        if drift_values:
            drift_data = [{'FIBER': fiber, 'DRIFT': drift} for fiber, drift in drift_values.items()]
            drift_table = Table(drift_data)
            ad[0].append(drift_table, extname='DRIFT')
            log.info(f'Added drift measurements to {ad.filename}')
        
        # Mark the file as having a dynamic wavelength solution
        ad.phu['DYNWVSOL'] = (True, 'Dynamic wavelength solution applied')
            

    
    return adinputs