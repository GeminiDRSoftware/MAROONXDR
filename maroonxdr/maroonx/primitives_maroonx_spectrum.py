"""
This module contains primitives to generate wavelength calibration solutions from reduced 1-D spectra.
"""
# ------------------------------------------------------------------------------
from astropy.table import Table
from recipe_system.utils.decorators import parameter_override
from . import parameters_maroonx_spectrum
from . import maroonx_etalon_fit
from . import maroonx_etalon_data
from .primitives_maroonx_echelle import MAROONXEchelle
from gempy.gemini import gemini_tools as gt
from geminidr.core import Spect
import multiprocessing
import numpy as np
import time
import traceback
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
    and utilized it to genterate mappings from pixel to wavelength
    (dynamic wavelength calibratiosn).
    """

    tagset = {"GEMINI", "MAROONX", "SPECT"}

    def __init__(self, adinputs, **kwargs):
        super(MaroonXSpectrum, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'maroonxdr.maroonx.lookups'
        self._param_update(parameters_maroonx_spectrum)

    def getPeaksAndPolynomials(self, adinputs=None, guess_file=None, 
                            fibers=[], orders=[], degree_sigma=4, degree_width=2, 
                            use_sigma_lr=True, show_plots=False,
                            plot_path="", multithreading=True, 
                            iterations = 3, **params):
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
            if fibers == []:
                log.warning("No fibers specified.  Finding all etalon fibers and extracting those.")
                # Read the HIERARCH FIBER keywords from the phu
                if ad.phu['HIERARCH FIBER1'] == 'Etalon':
                    fibers.append(1)
                if ad.phu['HIERARCH FIBER2'] == 'Etalon':
                    fibers.append(2)
                if ad.phu['HIERARCH FIBER3'] == 'Etalon':
                    fibers.append(3)
                if ad.phu['HIERARCH FIBER4'] == 'Etalon':
                    fibers.append(4)
                if ad.phu['HIERARCH FIBER5'] == 'Etalon':
                    fibers.append(5)

            if orders == []:
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
                for fiber, order, data, guess in maroonx_etalon_data.load_recordings(
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
                        maroonx_etalon_fit.iterative_fit,
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
                for fiber, order, data, guess in maroonx_etalon_data.load_recordings(
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
                    output = maroonx_etalon_fit.iterative_fit(
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
                peaks = maroonx_etalon_data.insert_peak_parameters(results)
                log.fullinfo(f"Adding poly extension to {ad.filename}")
                poly = maroonx_etalon_data.insert_polynomial_parameters(results)
                peaks_astropy = Table.from_pandas(peaks)
                poly_astropy = Table.from_pandas(poly)
                ad[0].PEAKS = peaks_astropy
                ad[0].POLY = poly_astropy

            if errors:
                log.warning("Errors were encountered during extraction")
                for fiber, order, msg in errors:
                    log.warning(f"Error extracting fiber {fiber} order {order}: {msg}")

        return adinputs
    