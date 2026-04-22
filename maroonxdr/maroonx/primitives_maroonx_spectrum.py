"""
This module contains primitives to generate wavelength
calibration solutions from reduced 1-D spectra.
"""

# ------------------------------------------------------------------------------
import multiprocessing
import os
import time
import traceback
import warnings
from copy import deepcopy
from pathlib import Path

import astrodata
import numpy as np
import pandas as pd
import scipy
from astropy import units
from astropy.table import Table
from astropy.time import Time, TimeDelta
from astroquery.simbad import Simbad
from geminidr.core import Spect
from gempy.adlibrary import dataselect
from gempy.gemini import gemini_tools as gt
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from maroonxdr.maroonx.maroonx_plots import (
    plot_residuals,
    plot_fiber_combination,
    plot_calibfiber_offset,
    plot_etalon_residuals,
    plot_exposuremeter,
)
from recipe_system.utils.decorators import parameter_override
from scipy.interpolate import interp1d
from scipy.signal import medfilt

from . import maroonx_utils, parameters_maroonx_spectrum
from .maroonx_echellespectrum.maroonxspectrum import MXSpectrum
from .maroonx_echellespectrum.wavelengthsolution import WavelengthSolution
from .maroonx_fit import maroonx_fit, set_logger
from .primitives_maroonx_echelle import MAROONXEchelle


# Monkey patch for barycorrpy compatibility with newer astroquery versions
def _patched_get_stellar_data(name=""):
    """
    Fixed version of barycorrpy.utils.get_stellar_data for newer astroquery.
    Based on solution from: https://github.com/shbhuk/barycorrpy/issues/59
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    warning = []

    customSimbad = Simbad()
    customSimbad.add_votable_fields(
        "ra", "dec", "pmra", "pmdec", "plx_value", "rvz_radvel"
    )

    obj = customSimbad.query_object(name)
    if obj is None:
        msg = (
            f"ERROR: {name} target not found. Check target name or "
            f"enter RA,Dec,PMRA,PMDec,Plx,RV,Epoch manually\n\n"
        )
        raise ValueError(msg)
    else:
        warning += [f"{name} queried from SIMBAD."]

    # Check for masked values
    if all([not x for x in [obj.mask[0][i] for i in obj.colnames]]) == False:
        warning += ["Masked values present in queried dataset"]

    obj = obj.filled(None)

    pos = SkyCoord(ra=obj["ra"], dec=obj["dec"], unit=(u.deg, u.deg))
    ra = pos.ra.value[0]
    dec = pos.dec.value[0]
    pmra = obj["pmra"][0]
    pmdec = obj["pmdec"][0]
    plx = obj["plx_value"][0]
    rv = obj["rvz_radvel"][0] * 1000  # SIMBAD output is in km/s. Converting to m/s
    epoch = 2451545.0

    star = {
        "ra": ra,
        "dec": dec,
        "pmra": pmra,
        "pmdec": pmdec,
        "px": plx,
        "rv": rv,
        "epoch": epoch,
    }

    # Fill Masked values with None. Again.
    for i in star:
        if star[i] > 1e10:
            star[i] = None

    warning += [f"Values queried from SIMBAD are {star}"]
    return star, warning


# Monkey patch the original function in barycorrpy
import barycorrpy.barycorrpy as bcp
import barycorrpy.utils as bcu

# update source module
bcu.get_stellar_data = _patched_get_stellar_data
# update alias in barycorrpy.barycorrpy
bcp.get_stellar_data = _patched_get_stellar_data

get_BC_vel = bcp.get_BC_vel
exposure_meter_BC_vel = bcp.exposure_meter_BC_vel

# ------------------------------------------------------------------------------


def _get_calibration_wavecal_path():
    """
    Get the path for the calibration flat file.
    Should probably be deprecated when dragons calib is implemented.
    """
    cwd = Path(os.getcwd())
    return cwd / 'calibrations' / 'processed_wavecal'


def _get_calibration_wavecal(adinputs):
    """
    Match and return calibration etalon file as astrodata object.
    Should probably be deprecated when dragons calib is implemented.
    """
    # Get the calibration etalons
    calib_path = _get_calibration_wavecal_path()

    arm_tag = "BLUE" if "BLUE" in adinputs[0].tags else "RED"
    etalons = dataselect.select_data(
        list(calib_path.glob("*.fits")), tags=[arm_tag, "ETALON", "PREPARED"]
    )
    ad_etalons = [astrodata.open(f) for f in etalons]

    adoutputs = []
    for ad in adinputs:
        science_time = ad.ut_datetime()

        # Find the etalon with minimum time difference
        def time_diff(etalon):
            return abs((etalon.ut_datetime() - science_time).total_seconds())

        closest_etalon = min(ad_etalons, key=time_diff)

        adoutputs.append(closest_etalon)
    return adoutputs


class LogExceptions:
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
            raise


@parameter_override
class MaroonXSpectrum(MAROONXEchelle, Spect):
    """
    This class contains primitives to reduce MAROON-X 1-D spectra.
    Code in this class takes already produced 1-D reduced spectra
    and utilizes it to generate mappings from pixel to wavelength
    (dynamic wavelength calibratiosn).
    """

    tagset = {"GEMINI", "MAROONX", "SPECT"}

    def _initialize(self, adinputs, **kwargs):
        super()._initialize(adinputs, **kwargs)
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
        fibers = params.get("fibers")
        if fibers is None:
            fibers = [1, 2, 3, 4, 5]

        for ad in adinputs:
            # Load static wavelength solution from the config
            statwavelength_file = maroonx_utils.get_statwavelength_filename(ad)
            log.info("Loading static wavelength file: %s", statwavelength_file)

            for fiber in range(1, 6):
                # Set up an initial value for all fibers
                setattr(ad[0], f"WLS_STATIC_FIBER_{fiber}", np.zeros((1, 1)))

            for fiber in fibers:
                # Check if the fiber is present in the data
                orders = getattr(ad[0], f"REDUCED_ORDERS_FIBER_{fiber}")
                if orders.size == 1:
                    # Save an empty array for this fiber
                    wls_data = np.zeros((1, 1))
                else:
                    # Load the static wavelength solution for this fiber
                    wls_dict = maroonx_utils.load_statwls_from_fits(
                        statwavelength_file, ext_name=f"FIBER_{fiber}", orders=orders
                    )
                    wls_data = np.vstack(
                        [wls_dict[str(int(order))] for order in orders]
                    )

                setattr(ad[0], f"WLS_STATIC_FIBER_{fiber}", wls_data)

            ad.update_filename(suffix=params["suffix"], strip=True)

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
        guess_file = params.get("guess_file")
        fibers = params.get("fibers")
        orders = params.get("orders")
        degree_sigma = params.get("degree_sigma")
        degree_width = params.get("degree_width")
        use_sigma_lr = params.get("use_sigma_lr")
        show_plots = params.get("show_plots", False)
        plot_path = params.get("plot_path", "")
        multithreading = params.get("multithreading")
        iterations = params.get("iterations")

        start_time = time.time()

        if len(adinputs) == 0:
            msg = "No input files"
            log.debug(msg)
            raise ValueError(msg)
        elif len(adinputs) >= 1:
            log.fullinfo(f"Extracting etalon lines from {len(adinputs)} files")

        for ad in adinputs:
            log.fullinfo(f"Extracting peaks from {ad.filename}")

            # See which fibers and orders we are extracting
            if fibers is None:
                log.warning(
                    "No fibers specified.  Finding all Etalon and LFC "
                    "fibers and extracting those."
                )
                fibers_list = []
                for fiber_num, fiber_type in enumerate(ad.fiber_setup(), start=1):
                    if fiber_type in ["Etalon", "LFC"]:
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
            set_logger(log)

            # Create worker pools for multithreading
            if multithreading:
                pool = multiprocessing.get_context("fork").Pool()
                log.fullinfo(f"Using {pool._processes} processes for extraction")
                for fiber, order, data, guess in maroonx_utils.load_recordings(
                    ad, guess_file, fibers, orders
                ):
                    log.fullinfo(f"{ad.filename} - Fitting fiber {fiber} order {order}")

                    # Remove pixels that are known to be bad
                    ############################
                    if order == 122:
                        # Remnant from the old pipeline, we should move
                        # this to the BPM at some point
                        data[1943] = np.nan
                        log.warning("Removed pixel 1943 in order 122")

                    if fiber == 5 and order == 94 and len(data) > 4000:
                        data[0:399] = 0
                        log.warning(
                            "Removed first 400 pixels in truncated "
                            "order 94 of fiber 5"
                        )
                    ############################

                    # Define callback functions to save results and errors
                    # for multithreading
                    def save_result(x, fiber=fiber, order=order):
                        x.fiber = fiber
                        x.order = order
                        x.recording_time = 0.0  # TODO extract from data set
                        results.append(x)

                    def error_callback(exc, fiber=fiber, order=order):
                        errors.append((fiber, order, str(exc)))
                        log.warning(
                            "%s - Error extracting fiber %s order %s: %s",
                            ad.filename, fiber, order, exc
                        )

                    if not (np.nanmedian(data) > 0.005):
                        log.warning(
                            "Skipped order %s for fiber %s for insufficient flux",
                            order, fiber
                        )
                        continue

                    # Asynchronously run iterative fit using the data
                    # yielded by the generator function
                    pool.apply_async(
                        maroonx_fit.iterative_fit,
                        kwds=dict(
                            input_spectrum=data,
                            degree_sigma=degree_sigma,
                            degree_width=degree_width,
                            iterations=iterations,
                            guess_spectrum=guess,
                            fiber=f"{fiber}_{order}",
                            plot_path=plot_path,
                            use_sigma_lr=use_sigma_lr,
                            show_plots=show_plots,
                        ),
                        callback=save_result,
                        error_callback=error_callback,
                    )

                pool.close()
                pool.join()

            # If multithreading is off, run in serial
            else:
                log.fullinfo("Using single process for extraction")
                for fiber, order, data, guess in maroonx_utils.load_recordings(
                    ad, guess_file, fibers, orders
                ):
                    log.fullinfo(f"{ad.filename} - Fitting fiber {fiber} order {order}")

                    # Remove pixels that are known to be bad
                    ############################
                    if order == 122:
                        # Remnant from the old pipeline, we should move
                        # this to the BPM at some point
                        data[1943] = np.nan
                        log.warning("Removed pixel 1943 in order 122")

                    if fiber == 5 and order == 94 and len(data) > 4000:
                        data[0:399] = 0
                        log.warning(
                            "Removed first 400 pixels in truncated "
                            "order 94 of fiber 5"
                        )
                    ############################

                    if not (np.nanmedian(data) > 0.005):
                        log.warning(
                            "Skipped order %s for fiber %s for insufficient flux",
                            order, fiber
                        )
                        continue

                    # Run iterative fit using the data yielded by the
                    # generator function in serial
                    try:
                        output = maroonx_fit.iterative_fit(
                            input_spectrum=data,
                            degree_sigma=degree_sigma,
                            degree_width=degree_width,
                            iterations=iterations,
                            guess_spectrum=guess,
                            fiber=f"{fiber}_{order}",
                            plot_path=plot_path,
                            use_sigma_lr=use_sigma_lr,
                            show_plots=show_plots,
                        )

                        # Save the results in a manner anologoous to the
                        # save_results callback function
                        output.fiber = fiber
                        output.order = order
                        output.recording_time = 0.0  # TODO extract from data set
                        results.append(output)
                    except Exception as exc:
                        errors.append((fiber, order, str(exc)))
                        log.warning(
                            "%s - Error extracting fiber %s order %s: %s",
                            ad.filename, fiber, order, exc
                        )

            # Record the time taken
            end_time = time.time()
            log.fullinfo(
                f"Finished extracting etalon lines in "
                f"{end_time - start_time:.2f} seconds"
            )

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
                    log.warning("Error extracting fiber %s order %s: %s", fiber, order, msg)

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
            Modified AstroData objects with new extensions and header keywords.

            **Extensions added:**

            - WLS_DYNAMIC_FIBER_* (1-5): 2D arrays of wavelength values (nm)
              indexed by [order, pixel] for each fiber
            - PEAK_DATA: Updated table with wavelength assignments and peak
              order numbers

            **Header keywords added:**

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
        fibers = params.get("fibers")
        symmetric_linefits = params.get("symmetric_linefits")
        n_knots = params.get("n_knots")
        thar = params.get("thar")
        ref_file = params.get("ref_file")
        report = params.get("report")


        for ad in adinputs:
            if "ETALON" not in ad.tags:
                log.warning(
                    "File %s is not ETALON. Skipping dynamic wavelength solution fitting.",
                    ad.filename
                )
                continue
 
            if report:
                # Create pdf for plots
                report_prefix = "spline_symmetrical_" if symmetric_linefits else "spline_"
                pdf = PdfPages(report_prefix + ad.filename.replace('.fits', '.pdf'))

            # Load the etalon spectrum
            mx_spectrum = MXSpectrum(ad, etalon_peaks_symmetric=symmetric_linefits)
            log.fullinfo(f"Processing etalon file: {ad.filename}")

            # Determine fibers to process if not provided
            if fibers is None:
                fibers = []
                for fiber_num, fiber_type in enumerate(ad.fiber_setup(), start=1):
                    if fiber_type == "Etalon":
                        fibers.append(fiber_num)
                log.fullinfo(f"Found etalon fibers: {fibers}")

            if ref_file is not None:
                # If a reference file and reference fiber is given, offset the etalon
                # position on all fibers by the offset found in the reference fiber.
                
                # TODO: Ask Andreas what this is supposed to do because currently
                # we do not know how this works in the old pipeline
                msg = "Reference file not implemented yet"
                log.debug(msg)
                raise NotImplementedError(msg)

            # Load reference wavelength solution from the config
            refwavelength_file = maroonx_utils.get_refwavelength_filename(ad)
            log.fullinfo(f"Loading reference wavelength file: {refwavelength_file}")

            # If chosen, apply ThAr wls to etalon frame
            if thar == True:
                ref_wavelength = {
                    1: maroonx_utils.load_refwls_from_fits(
                        refwavelength_file, ext_name="FIBER_2"
                    ),
                    2: maroonx_utils.load_refwls_from_fits(
                        refwavelength_file, ext_name="FIBER_2"
                    ),
                    3: maroonx_utils.load_refwls_from_fits(
                        refwavelength_file, ext_name="FIBER_3"
                    ),
                    4: maroonx_utils.load_refwls_from_fits(
                        refwavelength_file, ext_name="FIBER_4"
                    ),
                    5: maroonx_utils.load_refwls_from_fits(
                        refwavelength_file, ext_name="FIBER_4"
                    ),
                }
                for fiber in fibers:
                    wls_solution = WavelengthSolution(**ref_wavelength[fiber])
                    mx_spectrum.spectra[fiber].apply_wavelength_solution(wls_solution)
            else:
                for fiber in fibers:
                    log.fullinfo(
                        f"Apply static wavelength solution for peak number "
                        f"identification in fiber {fiber}."
                    )
                    mx_spectrum.spectra[fiber].apply_wavelength_vector()

            log.fullinfo(f"Apply etalon parameters from file {refwavelength_file}.")
            parameters = maroonx_utils.load_params_from_fits(
                refwavelength_file, ext_name="PARAMETERS"
            )
            for fiber in fibers:
                mx_spectrum.spectra[fiber].etalon_parameters = parameters

            # ========================================================= OK line

            # Guess the order #s of the etalon peak positions in the measured spectrum
            wl = None
            m = None
            mf = None
            o = None
            x = None

            inst_drift = None

            drifts = {}

            for fiber in fibers:
                log.info("Guess Etalon peak numbers for fiber %s", fiber)
                peak_data = mx_spectrum.spectra[fiber].guess_peak_numbers(debug=0)

                if report:
                    fig = mx_spectrum.spectra[fiber].plot_etalon_dispersion(
                        plot_title=f'(Fiber {fiber})', plot_mfraction=False)
                    pdf.savefig(fig)
                    plt.close('all')

                if wl is not None:
                    wl = np.concatenate((wl, peak_data["WAVELENGTH_BY_THAR"].values))
                    m = np.concatenate((m, peak_data["M"].values))
                    mf = np.concatenate((mf, peak_data["M_FRACTION"].values))
                    o = np.concatenate((o, peak_data["ORDER"].values))
                    x = np.concatenate((x, peak_data["CENTER"].values))
                else:
                    wl = peak_data["WAVELENGTH_BY_THAR"].values
                    m = peak_data["M"].values
                    mf = peak_data["M_FRACTION"].values
                    o = peak_data["ORDER"].values
                    x = peak_data["CENTER"].values

                # Calculate drift for this fiber
                residuals = (
                    _fc2min(
                        parameters,
                        peak_data["M"].values,
                        peak_data["WAVELENGTH_BY_THAR"].values,
                    )
                    / peak_data["WAVELENGTH_BY_THAR"].values
                    * 3e8
                )
                bad = np.where(
                    np.abs(residuals - np.nanmedian(residuals))
                    > 4.0 * np.nanstd(residuals)
                )
                residuals[bad] = np.nan
                if fiber == 5:
                    inst_drift = np.nanmean(residuals)
                    drifts[fiber] = inst_drift
                else:
                    drift = np.nanmean(residuals)
                    drifts[fiber] = drift

            # ==============================================================
            # This should probably be splitted into a separate primitive

            # Create new wls from etalon peaks based on the fitted
            # etalon gap size and dispersion model.
            wave = {}

            new_peak_data = []

            for fiber in fibers:
                log.info("Spline fit for fiber %s", fiber)
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
                        x = center.values
                        y = _peak_to_wavelength_spline(
                            speactra_peaks["M"][o], spectra.etalon_parameters
                        ).values
                    else:
                        x = (center.values)[::-1]
                        y = (
                            _peak_to_wavelength_spline(
                                speactra_peaks["M"][o], spectra.etalon_parameters
                            ).values
                        )[::-1]
                    knots = np.linspace(np.min(x) + 1, np.max(x) - 1, n_knots)
                    lsq = scipy.interpolate.LSQUnivariateSpline(x, y, knots, k=3)
                    r = (y - lsq(x)) / y * 3e8
                    good = np.where(np.abs(r) < 3.5 * np.nanstd(r))
                    lsq = scipy.interpolate.LSQUnivariateSpline(
                        x[good], y[good], knots, k=3, ext=3
                    )

                    xs_all = np.append(xs_all, x[good], axis=0)
                    orders_all = np.append(
                        orders_all, np.ones_like(x[good]) * o, axis=0
                    )
                    wavelengths_all = np.append(wavelengths_all, lsq(x[good]), axis=0)
                    residuals_all = np.append(
                        residuals_all, y[good] - lsq(x[good]), axis=0
                    )

                    xx = np.arange(len(spectra.data.wavelength.iloc[0]))
                    wavelengths = lsq(xx)

                    f = scipy.interpolate.interp1d(
                        x[good][-2:], y[good][-2:], fill_value="extrapolate"
                    )
                    wavelengths[wavelengths == np.max(wavelengths)] = f(
                        xx[wavelengths == np.max(wavelengths)]
                    )
                    f = scipy.interpolate.interp1d(
                        x[good][:2], y[good][:2], fill_value="extrapolate"
                    )
                    wavelengths[wavelengths == np.min(wavelengths)] = f(
                        xx[wavelengths == np.min(wavelengths)]
                    )

                    if fiber in wave:
                        wave[fiber].update({str(int(o)): wavelengths})
                    else:
                        wave[fiber] = {str(int(o)): wavelengths}

                # Plot residuals if report is True
                if report:
                    normalize_x = lambda x: x / np.max(x) * 2. - 1.
                    normalize_order = lambda o: (o - np.min(o)) / (np.max(o) - np.min(o)) * 2. - 1.
                    
                    fig = plot_residuals(
                        residuals_all, wavelengths_all,
                        normalize_order(orders_all), normalize_x(xs_all),
                        plottitle=f'Etalon residuals after Etalon-based spline fit (n_knots: {n_knots}) for fiber {fiber}',
                    )
                    pdf.savefig(fig)
                    plt.close('all')

                new_peak_data.append(speactra_peaks.copy())

                # Concatenate wavelengths arrays for each fiber and save them
                wave_arrays = [wave[fiber][o] for o in wave[fiber].keys()]
                setattr(ad[0], f"WLS_DYNAMIC_FIBER_{fiber}", np.stack(wave_arrays))

            for fiber in {1, 2, 3, 4, 5} - set(fibers):
                # If fiber is not in the list of fibers, save an empty array
                setattr(ad[0], f"WLS_DYNAMIC_FIBER_{fiber}", np.zeros((1, 1)))

            # Collect and reformat updated peak data
            new_peak_data = pd.concat(new_peak_data).set_index(
                ["FIBER", "ORDER", "M"], drop=False
            )
            new_peak_data.sort_index(level=[0, 1, 2], inplace=True)
            ad[0].PEAK_DATA = Table.from_pandas(new_peak_data)

            # Save meassured drift in header entries
            for fiber in fibers:
                ad[0].hdr[f"DRIFT_FIBER_{fiber}"] = (
                    round(drifts[fiber], 2),
                    "Drift in m/s",
                )
                log.info("Drift for fiber %s: %.2f m/s", fiber, drifts[fiber])

            ad.update_filename(suffix=params["suffix"], strip=True)
            if report:
                pdf.close()

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
        report : bool, optional
            If True, generate diagnostic plots of pixel shifts and spline fits.
        suffix : str, optional
            Suffix to append to output filenames. Default is empty string.

        Returns
        -------
        list of AstroData
            Modified AstroData objects with new extensions and header keywords.

            **Extensions added:**

            - WLS_SIMULTANEOUS_FIBER_* (fibers + ref_fiber): 2D arrays of
              drift-corrected wavelength values (nm) indexed by [order, pixel]

            **Header keywords added:**

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
        fibers = params.get("fibers", [2, 3, 4])
        ref_fiber = params.get("ref_fiber", 5)
        symmetric_linefits = params.get("symmetric_linefits", False)
        n_knots = params.get("n_knots", 30)
        wavecal = params.get("wavecal")
        report = params.get("report")

        # Resolve wavecal: use parameter or fall back to caldb
        if wavecal is None:
            wavecal_list = self.caldb.get_processed_wavecal(adinputs)
        else:
            wavecal_list = (wavecal, None)

        for science_ad, etalon_ad, _ in zip(
            *gt.make_lists(adinputs, *wavecal_list, force_ad=(1,))
        ):
            if etalon_ad is None:
                raise RuntimeError(
                    f"No processed wavecal listed for {science_ad.filename}"
                )
            log.fullinfo(f"Processing: {science_ad.filename} , {etalon_ad.filename}")
            log.fullinfo(f"Etalon reference fiber: {ref_fiber}")

            if report:
                # Create pdf for plots
                pdf_filename = science_ad.filename.replace('.fits', '_spline.pdf')
                pdf = PdfPages(pdf_filename)
                fig = None

            # Load the science and etalon spectrum
            science = MXSpectrum(science_ad, etalon_peaks_symmetric=symmetric_linefits)
            etalon = MXSpectrum(etalon_ad, etalon_peaks_symmetric=symmetric_linefits, wave_ext='WLS_DYNAMIC')

            # Load reference wavelength solution from the config
            refwavelength_file = maroonx_utils.get_refwavelength_filename(science_ad)
            parameters = maroonx_utils.load_params_from_fits(
                refwavelength_file, ext_name="PARAMETERS"
            )
            log.debug("Apply etalon parameters from file %s.", refwavelength_file)
            for fiber in fibers:
                etalon.spectra[fiber].etalon_parameters = parameters
            etalon.spectra[ref_fiber].etalon_parameters = parameters
            science.spectra[ref_fiber].etalon_parameters = parameters

            # Apply wavelength vector in etalon spectra to the etalon peak positions
            for fiber in fibers:
                etalon.spectra[fiber].apply_wavelength_vector()
            etalon.spectra[ref_fiber].apply_wavelength_vector()

            # # DEBUG TEST 1 output
            # debug_fibers = [2, 3, 4, 5]
            # debug_science = {}
            # debug_etalon = {}
            # for fiber in debug_fibers:
            #     sci_spec = science.spectra.get(fiber)
            #     eta_spec = etalon.spectra.get(fiber)
            #     if fiber in science.spectra and hasattr(science.spectra[fiber], 'peak_data'):
            #         debug_science[fiber] = {
            #             "peak_data": sci_spec.peak_data.to_dict(),
            #         }
            #     if fiber in etalon.spectra and hasattr(etalon.spectra[fiber], 'peak_data'):
            #         debug_etalon[fiber] = {
            #             "peak_data": eta_spec.peak_data.to_dict(),
            #         }
            # np.save(f"debug_test1_{science_ad.filename.replace('.fits', '.npy')}", debug_science, allow_pickle=True)
            # np.save(f"debug_test1_{etalon_ad.filename.replace('.fits', '.npy')}", debug_etalon, allow_pickle=True)
            # log.debug("DEBUG TEST 1: Saved debug_test1_science.npy and debug_test1_etalon.npy")
            # --------------


            # Calculate offsets in pixel space between science and etalon
            # frame for reference fiber
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
                    shift = x.loc[:] - x_ref.loc[:].reindex(
                        x.loc[:].index, method="nearest", tolerance=0.5
                    )
                    shift[
                        np.abs(shift)
                        > (np.abs(np.nanmedian(shift)) + 3 * np.nanstd(shift))
                    ] = np.nan
                    log.fullinfo(
                        f"Removed {np.count_nonzero(np.isnan(shift)):.0f} peaks "
                        f"in offset calculation for fiber {ref_fiber} in order {o}"
                    )
                    log.fullinfo(
                        f"Found {np.nanmean(shift):.3f} pix offset for fiber "
                        f"{ref_fiber} in order {o}"
                    )
                except KeyError:
                    log.warning(
                        "Could not reference etalon line lists for fiber %s in order %s",
                        ref_fiber, o
                    )

                if o == 94 and "RED" in etalon_ad.tags:
                    shift[0:600] = np.nanmedian(shift[600:])
                mask = np.isnan(shift)
                spl = scipy.interpolate.LSQUnivariateSpline(
                    shift.index[~mask], shift.values[~mask], [1000, 2000, 3000], k=3
                )

                shifts = np.append(
                    shifts,
                    shift.values
                    * etalon.spectra[ref_fiber].peak_data["DISPERSION_MPS"][o],
                )
                wavelengths = np.append(
                    wavelengths,
                    etalon.spectra[ref_fiber].peak_data["WAVELENGTH_BY_THAR"][o],
                )
                orders = np.append(orders, np.ones_like(shift.values, dtype=int) * o)
                x_refs = np.append(x_refs, x_ref.values)
                xs = np.append(xs, x.values)
                splfits = np.append(
                    splfits,
                    spl(shift.index)
                    * etalon.spectra[ref_fiber].peak_data["DISPERSION_MPS"][o],
                )

                # Apply offsets (smoothed with spline) to etalon frame
                # for science fibers (spline is evaluated at peak indices, not pixel positions)
                for fiber in fibers:
                    center_series = etalon.spectra[fiber].peak_data.loc[o, 'CENTER']
                    center_spl = spl(center_series.index)  # evaluate at pixel positions (index), not raw array

                    log.debug("BEFORE shift fiber %s order %s: %s", fiber, o, center_series.values[:5])
                    log.debug("APPLYING shift fiber %s order %s: %s", fiber, o, center_spl[:5])
                    etalon.spectra[fiber].peak_data.loc[o, 'CENTER'] = center_series.values - center_spl
                    log.debug("AFTER shift fiber %s order %s: %s", fiber, o, etalon.spectra[fiber].peak_data.loc[o, 'CENTER'].values[:5])

            # DEBUG TEST 2 output
            # 
            # --------------

            if report:
                fig = plot_calibfiber_offset(xs, shifts, orders, wavelengths, splfits,
                    fig=fig, plottitle=f'to {etalon_ad.filename}')

                pdf.savefig(fig)
                plt.close(fig)


            # Guess the order #s of the etalon peak positions in the measured
            # spectrum of the science fibers in the etalon spectrum
            # for fiber in fibers + [ref_fiber]:
            #     etalon.spectra[fiber].apply_wavelength_vector()
            #     etalon.spectra[fiber].guess_peak_numbers(debug=0)
            

            fig=None
            for fiber in fibers:
                etalon.spectra[fiber].apply_wavelength_vector()
                etalon_peak_data = etalon.spectra[fiber].guess_peak_numbers(debug=0, drop_outliers=False)
                if report:
                    residuals = _fc2min(
                        parameters, 
                        etalon_peak_data["M"].values, 
                        etalon_peak_data["WAVELENGTH_BY_THAR"].values) / etalon_peak_data["WAVELENGTH_BY_THAR"].values * 3e8
                    # residuals = _fc2min(
                    #     parameters,
                    #     etalon_peak_data["M"].values,
                    #     np.zeros_like(etalon_peak_data["WAVELENGTH_BY_THAR"].values))

                    bad = np.where(np.abs(residuals-np.nanmedian(residuals)) > 4.0 * np.nanstd(residuals))
                    residuals[bad] = np.nan
                    fig = plot_etalon_residuals(
                        etalon_peak_data["WAVELENGTH_BY_THAR"].values,
                        residuals,
                        etalon_peak_data["ORDER"].values,
                        plottitle=f'after sim cal correction (Fiber {fiber})',
                        plotnumber=fiber-2,
                        fig=fig,
                    )
                    
                    fig_debug = _plot_debug(fiber, etalon_peak_data, parameters, residuals)
                    pdf_debug = PdfPages(f"debug_fiber{fiber}.pdf")
                    pdf_debug.savefig(fig_debug)
                    pdf_debug.close()
                    plt.close(fig_debug)

            etalon.spectra[ref_fiber].apply_wavelength_vector()

            # Guess the order #s of the etalon peak positions in the measured
            # spectrum for reference fiber in science spectrum
            science.spectra[ref_fiber].apply_wavelength_vector()
            peak_data = science.spectra[ref_fiber].guess_peak_numbers(debug=0)

            residuals = (
                _fc2min(
                    parameters,
                    peak_data["M"].values,
                    peak_data["WAVELENGTH_BY_THAR"].values,
                )
                / peak_data["WAVELENGTH_BY_THAR"].values
                * 3e8
            )
            bad = np.where(
                np.abs(residuals - np.nanmedian(residuals)) > 4.0 * np.nanstd(residuals)
            )
            residuals[bad] = np.nan
            inst_drift = np.nanmean(residuals)
            
            if report:
                fig = plot_etalon_residuals(
                    peak_data["WAVELENGTH_BY_THAR"].values,
                    residuals,
                    peak_data["ORDER"].values,
                    plottitle=f'(Fiber {ref_fiber})',
                    plotnumber=ref_fiber-2,
                    fig=fig,
                )
                pdf.savefig(fig)

            # Create new wls from etalon peaks based on the given etalon gap
            # size and dispersion model.
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
                        x = center.values
                        y = _peak_to_wavelength_spline(
                            etalon.spectra[fiber].peak_data["M"][o],
                            etalon.spectra[fiber].etalon_parameters,
                        ).values
                    else:
                        x = (center.values)[::-1]
                        y = (
                            _peak_to_wavelength_spline(
                                etalon.spectra[fiber].peak_data["M"][o],
                                etalon.spectra[fiber].etalon_parameters,
                            ).values
                        )[::-1]
                    knots = np.linspace(min(x) + 1, max(x) - 1, n_knots)
                    mask = np.isnan(x)
                    if np.count_nonzero(mask) > 0:
                        log.warning(
                            "Found %s NANs in order %s for fiber %s while building the spline wls",
                            np.count_nonzero(mask), o, fiber
                        )
                    lsq = scipy.interpolate.LSQUnivariateSpline(
                        x[~mask], y[~mask], knots, k=3
                    )
                    r = (y - lsq(x)) / y * 3e8
                    good = np.where(np.abs(r - np.nanmedian(r)) < 3.5 * np.nanstd(r))
                    lsq = scipy.interpolate.LSQUnivariateSpline(
                        x[good], y[good], knots, k=3, ext=3
                    )

                    xs_all = np.append(xs_all, x[good], axis=0)
                    orders_all = np.append(
                        orders_all, np.ones_like(x[good]) * o, axis=0
                    )
                    wavelengths_all = np.append(wavelengths_all, lsq(x[good]), axis=0)
                    residuals_all = np.append(
                        residuals_all, y[good] - lsq(x[good]), axis=0
                    )

                    xx = np.arange(len(etalon.spectra[fiber].data.wavelength[92]))
                    wavelengths = lsq(xx)

                    f = scipy.interpolate.interp1d(
                        x[good][-2:], y[good][-2:], fill_value="extrapolate"
                    )
                    wavelengths[wavelengths == np.max(wavelengths)] = f(
                        xx[wavelengths == np.max(wavelengths)]
                    )
                    f = scipy.interpolate.interp1d(
                        x[good][:2], y[good][:2], fill_value="extrapolate"
                    )
                    wavelengths[wavelengths == np.min(wavelengths)] = f(
                        xx[wavelengths == np.min(wavelengths)]
                    )

                    f = f"fiber_{fiber}"
                    if f in wave:
                        wave[f].update({str(o): wavelengths})
                    else:
                        wave[f] = {str(o): wavelengths}

                if report:
                    normalize_x = lambda x: x / np.max(x) * 2. - 1.
                    normalize_order = lambda o: (o - np.min(o)) / (
                            np.max(o) - np.min(o)) * 2. - 1.
                    fig = plot_residuals(
                        residuals_all, wavelengths_all,
                        normalize_order(orders_all), normalize_x(xs_all),
                        plottitle=f'Etalon residuals after spline fit (n_knots: {n_knots}) for fiber {fiber}',
                    )
                    pdf.savefig(fig)            

            # Now for the reference fiber in the science spectrum. This should
            # be refactored at some point
            wavelengths_all = []
            residuals_all = []
            orders_all = []
            xs_all = []

            ref_spectra = science.spectra[ref_fiber]

            for o in ref_spectra.peak_data["CENTER"].index.levels[0]:
                if (
                    ref_spectra.peak_data["CENTER"][o].values[0]
                    < ref_spectra.peak_data["CENTER"][o].values[10]
                ):
                    x = ref_spectra.peak_data["CENTER"][o].values
                    y = _peak_to_wavelength_spline(
                        ref_spectra.peak_data["M"][o], ref_spectra.etalon_parameters
                    ).values
                else:
                    x = (ref_spectra.peak_data["CENTER"][o].values)[::-1]
                    y = (
                        _peak_to_wavelength_spline(
                            ref_spectra.peak_data["M"][o], ref_spectra.etalon_parameters
                        ).values
                    )[::-1]
                knots = np.linspace(min(x) + 1, max(x) - 1, n_knots)
                lsq = scipy.interpolate.LSQUnivariateSpline(x, y, knots, k=3)
                r = (y - lsq(x)) / y * 3e8
                good = np.where(np.abs(r - np.nanmedian(r)) < 3.5 * np.nanstd(r))
                try:
                    lsq = scipy.interpolate.LSQUnivariateSpline(
                        x[good], y[good], knots, k=3, ext=3
                    )
                except:
                    log.warning(
                        "Spline fit failed for reference fiber, order %s - try again", o
                    )
                    good = np.where(np.abs(r - np.nanmedian(r)) < 5 * np.nanstd(r))
                    try:
                        lsq = scipy.interpolate.LSQUnivariateSpline(
                            x[good], y[good], knots, k=3, ext=3
                        )
                        log.fullinfo(
                            "Spline fit succeeded for reference fiber, order %s with higher outlier threshold",
                            o
                        )
                    except:
                        log.warning(
                            "Spline fit failed again for reference fiber, order %s", o
                        )

                xs_all = np.append(xs_all, x[good], axis=0)
                orders_all = np.append(orders_all, np.ones_like(x[good]) * o, axis=0)
                wavelengths_all = np.append(wavelengths_all, lsq(x[good]), axis=0)
                residuals_all = np.append(residuals_all, y[good] - lsq(x[good]), axis=0)

                xx = np.arange(len(etalon.spectra[fiber].data.wavelength[92]))
                wavelengths = lsq(xx)

                f = scipy.interpolate.interp1d(
                    x[good][-2:], y[good][-2:], fill_value="extrapolate"
                )
                wavelengths[wavelengths == np.max(wavelengths)] = f(
                    xx[wavelengths == np.max(wavelengths)]
                )
                f = scipy.interpolate.interp1d(
                    x[good][:2], y[good][:2], fill_value="extrapolate"
                )
                wavelengths[wavelengths == np.min(wavelengths)] = f(
                    xx[wavelengths == np.min(wavelengths)]
                )

                f = f"fiber_{ref_fiber}"
                if f in wave:
                    wave[f].update({str(o): wavelengths})
                else:
                    wave[f] = {str(o): wavelengths}

            if report:
                fig = plot_residuals(
                    residuals_all, wavelengths_all,
                    normalize_order(orders_all), normalize_x(xs_all),
                    plottitle=f'Etalon residuals after spline fit (n_knots: {n_knots}) for fiber {ref_fiber}',
                )
                pdf.savefig(fig)
                pdf.close()
                log.stdinfo(f"Saved report file to {pdf_filename}")

            # Calculate mean drift in reference fiber
            rel_drift = np.nanmean(shifts)

            # Concatenate wavelengths arrays for each fiber and save them
            all_fibers = fibers + [ref_fiber]
            for fiber in all_fibers:
                f = f"fiber_{fiber}"
                wave_arrays = [wave[f][o] for o in wave[f].keys()]
                setattr(
                    science_ad[0],
                    f"WLS_SIMULTANEOUS_FIBER_{fiber}",
                    np.stack(wave_arrays),
                )

            for fiber in {1, 2, 3, 4, 5} - set(all_fibers):
                # If fiber is not in the list of fibers, save an empty array
                setattr(
                    science_ad[0], f"WLS_SIMULTANEOUS_FIBER_{fiber}", np.zeros((1, 1))
                )

            # Save meassured drift in header entries
            if inst_drift is not None:
                science_ad[0].hdr["INSTRUME_DRIFT"] = (
                    round(inst_drift, 2),
                    "Drift in m/s",
                )
                science_ad[0].hdr["RELATIVE_DRIFT"] = (
                    round(rel_drift, 2),
                    "Drift in m/s",
                )
                log.fullinfo(f"Instrument Drift: {inst_drift:.1f} m/s")
                log.fullinfo(
                    f"Relative drift measured in Fiber {ref_fiber}: {rel_drift:.1f} m/s"
                )
            log.fullinfo(f"Updated wavelength vector in {science_ad.filename}")

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
            - REDUCED_ORDERS_FIBER_6 (or _7):
              List of reduced orders.

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
        report = params.get("report")

        if combine_fibers is None:
            combine_fibers = [2, 3, 4]

        for ad in adinputs:

            if report:
                # Create pdf for plots
                tag = "fiber7_symmetrical" if symmetric_linefits else "fiber6"
                pdf_filename = ad.filename.replace('.fits', f'_{tag}.pdf')
                pdf = PdfPages(pdf_filename)

            # Get fiber data
            fiber_data = {}
            fiber_errors = {}
            fiber_wls = {}
            for fib in combine_fibers:
                fiber_data[fib] = getattr(ad[0], f"OPTIMAL_REDUCED_FIBER_{fib}")
                fiber_errors[fib] = getattr(ad[0], f"OPTIMAL_REDUCED_VAR_{fib}")

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

                intensity1 = fiber_data[2][order_idx].copy()
                intensity2 = fiber_data[3][order_idx].copy()
                intensity3 = fiber_data[4][order_idx].copy()
                
                scale1 = np.nansum(intensity1 / np.nansum(intensity2))
                intensity1 = intensity1 / scale1
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
                if "BLUE" in ad.tags:
                    log.fullinfo("Masking additional pixels in blue arm at pixel 196")
                    intensity1[196] = np.nan
                    intensity2[196] = np.nan
                    intensity3[196] = np.nan

                if "RED" in ad.tags:
                    log.fullinfo(
                        "Masking additional pixels in red arm at pixels 1793-1794"
                    )
                    intensity1[1793:1794] = np.nan
                    intensity2[1793:1794] = np.nan
                    intensity3[1793:1794] = np.nan

                mask1 = np.isnan(intensity1)
                mask2 = np.isnan(intensity2)
                mask3 = np.isnan(intensity3)

                try:
                    m = np.zeros_like(wave1, dtype=bool)
                    m[np.unique(wave1, return_index=True, return_inverse=True)[1]] = (
                        True
                    )
                    mask1[~m] = True

                    f1 = interp1d(
                        wave1[~mask1],
                        intensity1[~mask1],
                        fill_value="extrapolate",
                        kind="slinear",
                    )
                    f1e = interp1d(
                        wave1[~mask1],
                        error1[~mask1],
                        fill_value="extrapolate",
                        kind="slinear",
                    )
                    intensity1_2 = f1(wave2)
                    error1_2 = np.abs(f1e(wave2))
                except:
                    msg = (
                        f"Interpolation of flux in science fiber 1 of order "
                        f"{order} failed"
                    )
                    log.debug(msg)
                    intensity1_2 = intensity1
                    error1_2 = np.abs(error1) * np.nan

                try:
                    m = np.zeros_like(wave3, dtype=bool)
                    m[np.unique(wave3, return_index=True, return_inverse=True)[1]] = (
                        True
                    )
                    mask3[~m] = True

                    f3 = interp1d(
                        wave3[~mask3],
                        intensity3[~mask3],
                        fill_value="extrapolate",
                        kind="slinear",
                    )
                    f3e = interp1d(
                        wave3[~mask3],
                        error3[~mask3],
                        fill_value="extrapolate",
                        kind="slinear",
                    )
                    intensity3_2 = f3(wave2)
                    error3_2 = np.abs(f3e(wave2))
                except:
                    msg = (
                        f"Interpolation of flux in science fiber 3 of order "
                        f"{order} failed"
                    )
                    log.debug(msg)
                    intensity3_2 = intensity3
                    error3_2 = np.abs(error3) * np.nan

                intensity1_2[mask1] = np.nan
                intensity3_2[mask3] = np.nan

                median_intensity = np.nanmedian(
                    [intensity1_2, intensity2, intensity3_2], axis=0
                )

                weights1_2 = 1.0 / error1_2
                weights2 = 1.0 / error2
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

                clip1 = np.where(
                    np.nan_to_num(np.abs(intensity1_2 - median_intensity))
                    > kappa_sigma * np.nan_to_num(np.sqrt(error1_2))
                )
                clip2 = np.where(
                    np.nan_to_num(np.abs(intensity2 - median_intensity))
                    > kappa_sigma * np.nan_to_num(np.sqrt(error2))
                )
                clip3 = np.where(
                    np.nan_to_num(np.abs(intensity3_2 - median_intensity))
                    > kappa_sigma * np.nan_to_num(np.sqrt(error3_2))
                )

                clip1_n = np.size(clip1)
                clip2_n = np.size(clip2)
                clip3_n = np.size(clip3)

                while max([clip1_n, clip2_n, clip3_n]) > max_clips and kappa_sigma < 10:
                    log.warning(
                        "Number of maximum clipped pixels exceeded in order %s", order
                    )
                    kappa_sigma = kappa_sigma + 0.5

                    clip1 = np.where(
                        np.nan_to_num(np.abs(intensity1_2 - median_intensity))
                        > kappa_sigma * np.nan_to_num(np.sqrt(error1_2))
                    )
                    clip2 = np.where(
                        np.nan_to_num(np.abs(intensity2 - median_intensity))
                        > kappa_sigma * np.nan_to_num(np.sqrt(error2))
                    )
                    clip3 = np.where(
                        np.nan_to_num(np.abs(intensity3_2 - median_intensity))
                        > kappa_sigma * np.nan_to_num(np.sqrt(error3_2))
                    )

                    clip1_n = np.size(clip1)
                    clip2_n = np.size(clip2)
                    clip3_n = np.size(clip3)

                log.fullinfo("Kappa_sigma in order %s: %.1f", order, kappa_sigma)
                log.fullinfo("Clipped %s pixels in fiber 2 of order %s", clip1_n, order)
                log.fullinfo("Clipped %s pixels in fiber 3 of order %s", clip2_n, order)
                log.fullinfo("Clipped %s pixels in fiber 4 of order %s", clip3_n, order)

                weights1_2[clip1] = 0
                weights2[clip2] = 0
                weights3_2[clip3] = 0

                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", r"All-NaN slice encountered")

                    weights = np.nansum([weights1_2, weights2, weights3_2], axis=0)
                    bad = np.where(weights == 0)
                    weights1_2[bad] = np.nan
                    weights2[bad] = np.nan
                    weights3_2[bad] = np.nan

                    intensity = np.average(
                        [
                            np.nan_to_num(intensity1_2),
                            np.nan_to_num(intensity2),
                            np.nan_to_num(intensity3_2),
                        ],
                        weights=[weights1_2, weights2, weights3_2],
                        axis=0,
                    )
                    error = np.abs(
                        1.0 / np.sum([weights1_2, weights2, weights3_2], axis=0)
                    )

                mask = np.isnan(intensity)
                if 1000 > np.count_nonzero(mask) >= 0:
                    log.debug(
                        "Number of NANs fixed in order %s: %s",
                        order, np.count_nonzero(mask)
                    )
                    f = interp1d(
                        wave2[~mask],
                        intensity[~mask],
                        fill_value="extrapolate",
                        kind="slinear",
                    )
                    intensity = f(wave2)
                else:
                    log.warning(
                        "Too many NANs found in order %s: %s",
                        order, np.count_nonzero(mask)
                    )

                error[mask] = 1e6
                error[np.isnan(error)] = 1e6

                combined_wavelength.update({order: wave2})
                combined_intensity.update({order: intensity})
                combined_error.update({order: error})

                # Generate diagnostic plot for this order
                if report:
                    fig = plot_fiber_combination(
                        order, wave1, wave2, wave3,
                        intensity1_2, intensity2, intensity3_2,
                        weights1_2, weights2, weights3_2,
                        intensity, error,
                        error1_2, error2, error3_2,
                        median_intensity, kappa_sigma
                    )
                    pdf.savefig(fig)
                    plt.close(fig)

            # Close PDF after all orders
            if report:
                pdf.close()
                log.stdinfo(f"Saved diagnostic report to {pdf_filename}")

            # ====================================================================

            # Store combined results as new extensions
            combined_intensity_array = np.vstack(
                [combined_intensity[o] for o in orders]
            )
            combined_error_array = np.vstack([combined_error[o] for o in orders])
            combined_wavelength_array = np.vstack(
                [combined_wavelength[o] for o in orders]
            )

            # Define final combined fiber number
            target_fiber = 6 if not symmetric_linefits else 7

            setattr(
                ad[0], f"OPTIMAL_REDUCED_FIBER_{target_fiber}", combined_intensity_array
            )
            setattr(ad[0], f"OPTIMAL_REDUCED_VAR_{target_fiber}", combined_error_array)
            setattr(ad[0], f"REDUCED_ORDERS_FIBER_{target_fiber}", orders)

            if symmetric_linefits:
                setattr(
                    ad[0],
                    f"WLS_SIMULTANEOUS_SYM_FIBER_{target_fiber}",
                    combined_wavelength_array,
                )
            else:
                setattr(
                    ad[0],
                    f"WLS_SIMULTANEOUS_FIBER_{target_fiber}",
                    combined_wavelength_array,
                )

            # Mark history
            fiber_list = ",".join(map(str, combine_fibers))
            gt.mark_history(
                ad,
                primname=self.myself(),
                keyword="FIBER_COMBINATION",
                comment=f"combined_fibers_{fiber_list}_to_{target_fiber}",
            )

            log.fullinfo(
                f"Combined fibers {fiber_list} into fiber {target_fiber} "
                f"for {ad.filename}"
            )
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
            exposure meter statistics added to the first extension header.

            **BERV values (m/s):**

            - BERV_MIDPOINT: BERV at nominal exposure midpoint
            - BERV_FLUXWEIGHTED_PC: BERV at flux-weighted midpoint (PC channel)
            - BERV_FLUXWEIGHTED_FRD: BERV at flux-weighted midpoint (FRD channel)
            - BERV_DIFFERENCE_PC: Difference between flux-weighted and nominal
            - BERV_DIFFERENCE_FRD: Difference between flux-weighted and nominal

            **Timing information:**

            - UTC_START: Corrected UTC start time (ISO format)
            - UTC_MIDPOINT: UTC midpoint time
            - UTC_FLUXWEIGHTED_PC: Flux-weighted UTC (PC channel)
            - UTC_FLUXWEIGHTED_FRD: Flux-weighted UTC (FRD channel)
            - UTC_CORRECTION: Applied time correction in seconds
            - JD_UTC_START: Julian date at start
            - JD_UTC_MIDPOINT: Julian date at midpoint
            - JD_UTC_FLUXWEIGHTED_PC: Flux-weighted JD (PC channel)
            - JD_UTC_FLUXWEIGHTED_FRD: Flux-weighted JD (FRD channel)

            **Exposure meter statistics (counts):**

            - COUNTS_PC_MIN/MAX/MEDIAN/STD: PC channel statistics
            - COUNTS_FRD_MIN/MAX/MEDIAN/STD: FRD channel statistics
            - COUNTS_PC_ZP: Applied PC zeropoint
            - COUNTS_FRD_ZP: Applied FRD zeropoint
            - SCALEFACTOR: Ratio of FRD to PC median counts

            **Target information:**

            - BERV_SIMBAD_TARGET: Target name used for BERV calculation

        References
        ----------
        .. [1] barycorrpy: https://github.com/shbhuk/barycorrpy
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Get parameters
        target_name = params.get("target_name")
        simbad_target_name = params.get("simbad_target_name")
        use_coords = params.get("use_coords")
        zp_pc = params.get("zp_pc")
        zp_frd = params.get("zp_frd")
        start_time = params.get("start_time")
        report = params.get("report")

        # Cache stellar data per target name to avoid repeated SIMBAD TAP hits.
        # Uses bcp.get_stellar_data by reference so it works with or without the
        # monkey-patch applied above.
        stellar_data_cache = {}

        for ad in adinputs:
            # Read target name from header
            target = ad.object()
            log.debug("%s, %s", target_name, target)
            # Handle user-supplied target name logic
            if target_name is not None:
                if target_name in target:
                    log.fullinfo("Selected target name: %s", target)
                    if simbad_target_name and len(simbad_target_name) > 0:
                        log.fullinfo(
                            "Replaced with user supplied name: %s", simbad_target_name
                        )
                        target = simbad_target_name
                else:
                    log.warning("Skip file %s", ad.filename)
                    continue

            # Fetch stellar data once per unique target name
            if not use_coords and target not in stellar_data_cache:
                try:
                    star, _ = bcp.get_stellar_data(target)
                    stellar_data_cache[target] = star
                    log.fullinfo("Fetched stellar data for %s from SIMBAD", target)
                except Exception as e:
                    log.warning("Target %s not recognized by SIMBAD: %s", target, e)
                    log.warning(
                        "Skip file %s - consider using --use_coords option", ad.filename
                    )
                    continue

            # Calculate time correction based on selected method
            utc_start = Time(ad.ut_datetime(), format="datetime", scale="utc")
            time_correction = TimeDelta(0, format="sec")

            mjd = ad[0].telescope_mjd(pretty=True)
            exptime = ad[0].exposure_time(pretty=True)
            if start_time == "mjd_start":
                time_correction = mjd - utc_start
            elif start_time == "mjd_end":
                if "BLUE" in ad.tags:
                    # 52 and 100 come from legacy code, no docs about them
                    offset = exptime + TimeDelta(52.0, format="sec")
                else:
                    offset = exptime + TimeDelta(100.0, format="sec")
                time_correction = (mjd - offset) - utc_start

            # Extract times and exposure meter readings for a given exposure
            emeter = self._exposuremeterStats(
                ad, utc_start, exptime, zp_pc=zp_pc, zp_frd=zp_frd
            )
            log.fullinfo(f'Exposuremeter stats (PC channel): {emeter["pc"]["stats"]}')
            log.fullinfo(f'Exposuremeter stats (FRD channel): {emeter["frd"]["stats"]}')

            # Calculate BERV at start, mid, end of exposure
            utc_start += time_correction
            utc_mid = utc_start + exptime / 2
            utc_end = utc_start + exptime

            # Call barycorrpy for barycentric calculations
            barycorrpy_kwargs = {
                "lat": 19.823801,
                "longi": -155.469047,
                "alt": 4213,
                "ephemeris": "de430",
                "zmeas": 0.0,
                "predictive": False,
                "leap_update": False,
            }
            if use_coords:
                ra = ad[0].hdr.get("MAROONX TELESCOPE TELRA") * 15.0
                dec = ad[0].hdr.get("MAROONX TELESCOPE TELDEC")
                barycorrpy_kwargs["ra"] = ra
                barycorrpy_kwargs["dec"] = dec
                log.fullinfo(f"Using RA: {ra:.4f} deg, DEC: {dec:.4f} deg")
            else:
                star = stellar_data_cache[target]
                barycorrpy_kwargs.update({
                    "ra": star["ra"],
                    "dec": star["dec"],
                    "pmra": star.get("pmra"),
                    "pmdec": star.get("pmdec"),
                    "px": star.get("px"),
                    "rv": star.get("rv"),
                    "epoch": star.get("epoch", 2451545.0),
                })
                log.fullinfo(f"Using target name: {target}")

            # BVC for nominal exposure midtime
            result1, warning, status = get_BC_vel(
                JDUTC=[utc_start.jd, utc_mid.jd, utc_end.jd], **barycorrpy_kwargs
            )

            # Average BVC over exposure weighted by exposure meter readings - BEST
            if emeter["pc"]["times"] is not None:
                result3, JDUTCMID_pc, warning3, status3 = exposure_meter_BC_vel(
                    JDUTC=emeter["pc"]["times"].jd + time_correction.jd,
                    expmeterflux=emeter["pc"]["readings"],
                    **barycorrpy_kwargs,
                )
            if emeter["frd"]["times"] is not None:
                result4, JDUTCMID_frd, warning4, status4 = exposure_meter_BC_vel(
                    JDUTC=emeter["frd"]["times"].jd + time_correction.jd,
                    expmeterflux=emeter["frd"]["readings"],
                    **barycorrpy_kwargs,
                )

            # dvdt values in m/s/s
            m_s = units.Unit("m/s")
            m_s_s = units.Unit("m/s/s")
            BC_dvdt = [
                (result1[2] - result1[0]) * m_s / (exptime),
                (result1[2] - result1[1]) * m_s / (exptime / 2),
                (result1[1] - result1[0]) * m_s / (exptime / 2),
            ]

            # Mid values
            if emeter["pc"]["times"] is not None:
                utc_fluxmid_pc = Time(JDUTCMID_pc, format="jd", scale="utc")
                utc_fluxmid_frd = Time(JDUTCMID_frd, format="jd", scale="utc")
                BC_fluxmid_pc = result3
                BC_fluxmid_frd = result4

            # If times are None, there was a gap in the photometer data and midpoints
            # and BERVs are taken from the nominal midpoint.
            if emeter["pc"]["times"] is None:
                utc_fluxmid_pc = utc_mid
                BC_fluxmid_pc = result1[1]

            if emeter["frd"]["times"] is None:
                utc_fluxmid_frd = utc_mid
                BC_fluxmid_frd = result1[1]

            # Log results
            log.fullinfo(f"File: {ad.filename}")
            log.fullinfo(f"Target: {target}")
            log.fullinfo(f"Exptime: {exptime} sec")
            log.fullinfo(
                f"BERV dv/dt: "
                f"{BC_dvdt[0].to(m_s_s).value:.4f}, "
                f"{BC_dvdt[1].to(m_s_s).value:.4f}, "
                f"{BC_dvdt[2].to(m_s_s).value:.4f} m/s/s"
            )
            log.fullinfo(f"BERV_MIDPOINT:           {result1[1]:.2f} m/s")
            log.fullinfo(f"BERV_FLUXWEIGHTED_PC:    {BC_fluxmid_pc:.2f} m/s")
            log.fullinfo(f"BERV_FLUXWEIGHTED_FRD:   {BC_fluxmid_frd:.2f} m/s")

            # Scale factor
            scale_factor = (
                emeter["frd"]["stats"]["median"] / emeter["pc"]["stats"]["median"]
            )

            # Save header entries
            ad[0].hdr.set("UTC_START", f"{utc_start.isot}", "UTC time at start of exposure")
            ad[0].hdr.set("UTC_CORRECTION", f"{time_correction.sec:.1f}", "UTC correction applied [s]")
            ad[0].hdr.set("UTC_MIDPOINT", f"{utc_mid.isot}", "UTC time at nominal midpoint of exposure")
            ad[0].hdr.set("UTC_FLUXWEIGHTED_PC", f"{utc_fluxmid_pc.isot}", "UTC flux-weighted midpoint PC")
            ad[0].hdr.set("UTC_FLUXWEIGHTED_FRD", f"{utc_fluxmid_frd.isot}", "UTC flux-weighted midpoint FRD")
            ad[0].hdr.set("JD_UTC_START", f"{utc_start.jd:.7f}", "JD (UTC) at start of exposure")
            ad[0].hdr.set("JD_UTC_MIDPOINT", f"{utc_mid.jd:.7f}", "JD (UTC) at nominal midpoint of exposure")
            ad[0].hdr.set("JD_UTC_FLUXWEIGHTED_PC", f"{utc_fluxmid_pc.jd:.7f}", "JD (UTC) flux-weighted midpoint PC")
            ad[0].hdr.set(
                "JD_UTC_FLUXWEIGHTED_FRD", f"{utc_fluxmid_frd.jd:.7f}", "JD (UTC) flux-weighted midpoint FRD"
            )

            ad[0].hdr.set("BERV_SIMBAD_TARGET", target, "Target name used for BERV lookup")
            ad[0].hdr.set("BERV_MIDPOINT", f"{result1[1]:.2f}", "BERV at nominal midpoint [m/s]")
            ad[0].hdr.set("BERV_FLUXWEIGHTED_PC", f"{BC_fluxmid_pc:.2f}", "BERV at flux-weighted midpoint PC [m/s]")
            ad[0].hdr.set("BERV_FLUXWEIGHTED_FRD", f"{BC_fluxmid_frd:.2f}", "BERV at flux-weighted midpoint FRD [m/s]")
            ad[0].hdr.set(
                "BERV_DIFFERENCE_PC", f"{(BC_fluxmid_pc - result1[1]):.2f}", "BERV difference flux-weighted minus midpoint PC [m/s]"
            )
            ad[0].hdr.set(
                "BERV_DIFFERENCE_FRD", f"{(BC_fluxmid_frd - result1[1]):.2f}", "BERV difference flux-weighted minus midpoint FRD [m/s]"
            )

            ad[0].hdr.set(
                "COUNTS_PC_MIN", f"{emeter['pc']['stats']['min']:.2f}", "PC exposure meter counts minimum"
            )
            ad[0].hdr.set(
                "COUNTS_PC_MAX", f"{emeter['pc']['stats']['max']:.2f}", "PC exposure meter counts maximum"
            )
            ad[0].hdr.set(
                "COUNTS_PC_MEDIAN", f"{emeter['pc']['stats']['median']:.2f}", "PC exposure meter counts median"
            )
            ad[0].hdr.set(
                "COUNTS_PC_STD", f"{emeter['pc']['stats']['std']:.2f}", "PC exposure meter counts std dev"
            )
            ad[0].hdr.set(
                "COUNTS_PC_ZP", f"{emeter['pc']['stats']['zeropoint']:.2f}", "PC exposure meter zero point"
            )

            ad[0].hdr.set(
                "COUNTS_FRD_MIN", f"{emeter['frd']['stats']['min']:.2f}", "FRD exposure meter counts minimum"
            )
            ad[0].hdr.set(
                "COUNTS_FRD_MAX", f"{emeter['frd']['stats']['max']:.2f}", "FRD exposure meter counts maximum"
            )
            ad[0].hdr.set(
                "COUNTS_FRD_MEDIAN", f"{emeter['frd']['stats']['median']:.2f}", "FRD exposure meter counts median"
            )
            ad[0].hdr.set(
                "COUNTS_FRD_STD", f"{emeter['frd']['stats']['std']:.2f}", "FRD exposure meter counts std dev"
            )
            ad[0].hdr.set(
                "COUNTS_FRD_ZP", f"{emeter['frd']['stats']['zeropoint']:.2f}", "FRD exposure meter zero point"
            )
            ad[0].hdr.set("SCALEFACTOR", f"{scale_factor:.1f}", "Flux scale factor FRD/PC")

            log.fullinfo(
                f"Barycentric velocity calculation completed for {ad.filename}"
            )

            if report:
                dt1 = TimeDelta(exptime.sec, format="sec")
                utc_end = utc_start + dt1
                fig = plot_exposuremeter(
                    emeter["pc"]["context"]["times"],
                    emeter["pc"]["context"]["readings"],
                    emeter["frd"]["context"]["times"],
                    emeter["frd"]["context"]["readings"],
                    utc_start,
                    utc_end,
                    emeter["pc"]["stats"]["zeropoint"],
                    emeter["frd"]["stats"]["zeropoint"],
                    target,
                    exptime.sec,
                )
                pdf_filename = ad.filename.replace('.fits', '_exposuremeter.pdf')
                with PdfPages(pdf_filename) as pdf:
                    pdf.savefig(fig)
                plt.close(fig)
                log.fullinfo("Saved exposure meter diagnostic plot to %s", pdf_filename)

            ad.update_filename(suffix=params["suffix"], strip=True)

        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        log.debug(gt.log_message("primitive", self.myself(), "complete"))
        return adinputs

    def _extract_spectra_data(self, ad, ext_index, fibers, show_wavelength):
        """
        Helper function to extract spectrum data from an AstroData object.

        Parameters
        ----------
        ad : AstroData
            AstroData object with extracted spectra
        ext_index : int
            Extension index to extract from (0 for Blue, 1 for Red in bundles)
        fibers : list of int
            Fiber numbers to extract
        show_wavelength : bool
            Whether to extract wavelength solutions

        Returns
        -------
        dict
            Dictionary mapping fiber numbers to spectrum data
        """
        log = self.log
        spectra_data = {}

        # Validate extension index
        if ext_index >= len(ad):
            log.warning(f"Extension {ext_index} not found in {ad.filename}")
            return spectra_data

        for fiber in fibers:
            # Get both optimal and box extractions if available
            optimal_flux = getattr(ad[ext_index], f"OPTIMAL_REDUCED_FIBER_{fiber}", None)
            box_flux = getattr(ad[ext_index], f"BOX_REDUCED_FIBER_{fiber}", None)

            # Filter out placeholder (1,1) arrays - these indicate fiber was not extracted
            if optimal_flux is not None and optimal_flux.shape == (1, 1):
                optimal_flux = None
            if box_flux is not None and box_flux.shape == (1, 1):
                box_flux = None

            if optimal_flux is None and box_flux is None:
                log.warning(
                    f"Fiber {fiber}: No valid extraction data found in {ad.filename}, skipping"
                )
                continue

            # Get order numbers
            orders = getattr(ad[ext_index], f"REDUCED_ORDERS_FIBER_{fiber}", None)
            if orders is None or len(orders) <= 1:
                log.warning(
                    f"Fiber {fiber}: No valid order data found in {ad.filename} ext {ext_index}, skipping"
                )
                continue

            # Get wavelength solution if requested
            wavelength = None
            if show_wavelength:
                # Try dynamic first, then static
                wavelength = getattr(ad[ext_index], f"WLS_SIMULTANEOUS_FIBER_{fiber}", None)
                if wavelength is None:
                    wavelength = getattr(ad[ext_index], f"WLS_STATIC_FIBER_{fiber}", None)

                if wavelength is None:
                    log.warning(
                        f"Fiber {fiber}: No wavelength solution found, "
                        f"displaying in pixel space"
                    )

            spectra_data[fiber] = {
                "optimal_flux": optimal_flux,
                "box_flux": box_flux,
                "orders": orders,
                "wavelength": wavelength,
            }

            # Log available extractions
            available_extractions = []
            if optimal_flux is not None:
                available_extractions.append("optimal")
            if box_flux is not None:
                available_extractions.append("box")

            log.fullinfo(
                f"Fiber {fiber}: {len(orders)} orders, "
                f"wavelength: {wavelength is not None}, "
                f"extractions: {', '.join(available_extractions)}"
            )

        return spectra_data

    def displaySpectra(self, adinputs=None, **params):
        """
        Display extracted spectra in browser using interactive Bokeh visualization.

        This primitive launches an interactive Bokeh server that displays MaroonX
        extracted spectra in a web browser. Users can zoom, pan, and inspect
        individual orders for quality assessment. The primitive blocks execution
        until the user clicks "Submit" or closes the browser window.

        The visualizer automatically opens at http://localhost:5006 (or next
        available port) and provides interactive plots showing all extracted
        orders for selected fibers. If wavelength calibration is available and
        show_wavelength=True, spectra are displayed in wavelength space;
        otherwise pixel space is used.

        Parameters
        ----------
        adinputs : list of AstroData
            Input AstroData objects with extracted spectra. Must have
            OPTIMAL_REDUCED_FIBER_* or BOX_REDUCED_FIBER_* extensions.
        fibers : list of int
            Fiber numbers to display (e.g., [2, 3, 4] for science fibers,
            [6] for combined fiber, [5] for calibration fiber).
            Default is [2, 3, 4].
        show_wavelength : bool
            If True and wavelength solution exists (WLS_STATIC_FIBER_* or
            WLS_DYNAMIC_FIBER_* extensions), display spectra vs wavelength (nm).
            If False or no wavelength solution, display vs pixel number.
            Default is False.

        Returns
        -------
        list of AstroData
            Unmodified input AstroData objects (this is a visualization-only
            primitive with no data modification).

        Notes
        -----
        The Bokeh server runs in a separate thread and communicates with the
        browser via websockets. Browser opens automatically based on user
        configuration in ~/.dragons/dragonsrc:

        [interactive]
        browser = chrome          # or firefox, safari
        theme = dark_minimal      # or light_minimal
        port_number = 5006        # starting port

        If the specified port is occupied, the server automatically tries
        subsequent ports up to 65535.

        The visualizer provides:
        - Dropdown to select fiber
        - Dropdown to select individual orders or "All Orders" view
        - Dropdown to select extraction type (optimal or box)
        - Interactive Bokeh plot with zoom, pan, hover tooltips
        - Submit button to continue reduction

        Examples
        --------
        Display science fibers after extraction (pixel space):
        >>> p.displaySpectra(fibers=[2, 3, 4])

        Display calibrated spectra in wavelength space:
        >>> p.displaySpectra(fibers=[2, 3, 4], show_wavelength=True)

        Display combined fiber spectrum:
        >>> p.displaySpectra(fibers=[6], show_wavelength=True)
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        # Get parameters
        fibers = params["fibers"]
        show_wavelength = params["show_wavelength"]

        if fibers is None:
            fibers = [2, 3, 4]

        log.fullinfo(f"Selected fibers: {fibers}")

        # Import Bokeh visualizer components
        try:
            from geminidr.interactive import server
            from .interactive.spectrum_viewer import SpectrumVisualizer
        except ImportError as e:
            log.warning(
                f"Failed to import Bokeh interactive components: {e}\n"
                f"Skipping spectrum display. Install bokeh to enable visualization."
            )
            return adinputs

        # Process all AstroData objects and collect spectra by file and arm
        all_files_data = []

        for ad in adinputs:
            log.fullinfo(f"Preparing spectrum display for {ad.filename}")

            # Determine if this is a bundle or single arm
            if 'BUNDLE' in ad.tags:
                # Bundle: ext 0 = Blue, ext 1 = Red
                log.fullinfo(f"{ad.filename} is a BUNDLE (ext 0=Blue, ext 1=Red)")

                arms_data = {}

                # Extract Blue arm (extension 0)
                blue_spectra = self._extract_spectra_data(ad, 0, fibers, show_wavelength)
                if blue_spectra:
                    arms_data['Blue'] = blue_spectra
                else:
                    log.warning(f"No valid Blue arm spectra in {ad.filename}")

                # Extract Red arm (extension 1)
                if len(ad) > 1:
                    red_spectra = self._extract_spectra_data(ad, 1, fibers, show_wavelength)
                    if red_spectra:
                        arms_data['Red'] = red_spectra
                    else:
                        log.warning(f"No valid Red arm spectra in {ad.filename}")
                else:
                    log.warning(f"{ad.filename} marked as BUNDLE but has only 1 extension")

                if arms_data:
                    all_files_data.append({
                        'filename': ad.filename,
                        'arms_data': arms_data
                    })

            elif 'BLUE' in ad.tags:
                # Single Blue arm
                log.fullinfo(f"{ad.filename} is a single Blue arm")
                blue_spectra = self._extract_spectra_data(ad, 0, fibers, show_wavelength)
                if blue_spectra:
                    all_files_data.append({
                        'filename': ad.filename,
                        'arms_data': {'Blue': blue_spectra}
                    })
                else:
                    log.warning(f"No valid spectra in {ad.filename}")

            elif 'RED' in ad.tags:
                # Single Red arm
                log.fullinfo(f"{ad.filename} is a single Red arm")
                red_spectra = self._extract_spectra_data(ad, 0, fibers, show_wavelength)
                if red_spectra:
                    all_files_data.append({
                        'filename': ad.filename,
                        'arms_data': {'Red': red_spectra}
                    })
                else:
                    log.warning(f"No valid spectra in {ad.filename}")

            else:
                # Unknown type - try to extract from first extension
                log.fullinfo(f"{ad.filename} has no BUNDLE/BLUE/RED tag, treating as single arm")
                spectra = self._extract_spectra_data(ad, 0, fibers, show_wavelength)
                if spectra:
                    all_files_data.append({
                        'filename': ad.filename,
                        'arms_data': {'Data': spectra}
                    })
                else:
                    log.warning(f"No valid spectra in {ad.filename}")

        if not all_files_data:
            log.warning("No valid spectra found in any input files, skipping display")
            return adinputs

        # Create visualizer for all files
        log.stdinfo(f"Launching interactive spectrum viewer for {len(all_files_data)} file(s)")
        log.stdinfo(f"Browser will open automatically at http://localhost:5006")
        log.stdinfo(f"Click 'Submit' or close browser to continue reduction")

        visualizer = SpectrumVisualizer(
            all_files_data,
            title="MaroonX Extracted Spectra",
            primitive_name="displaySpectra",
        )

        # Launch browser (blocks until user submits)
        try:
            server.interactive_fitter(visualizer)
            log.fullinfo("User submitted spectrum viewer, continuing reduction")
        except KeyboardInterrupt:
            log.warning("User aborted spectrum viewer")
            raise

        log.debug(gt.log_message("primitive", self.myself(), "complete"))
        return adinputs

    def separateArmStreams(self, adinputs=None, **params):
        """
        Separate input AstroData objects into blue and red arm streams.

        This primitive takes a list of AstroData objects and separates them
        into two lists based on their tags: one for the blue arm and one for
        the red arm.

        Parameters
        ----------
        adinputs : list of AstroData
            Input list of AstroData objects to be separated.

        Returns
        -------
        list of AstroData
            List of AstroData objects tagged as 'BLUE' arm. The RED arm
            objects are stored in self.streams['RED'].
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        # Initialize dictionary to beindexed by ARCHNAME and arm tag
        arm_dict = {}
        for ad in adinputs:
            archname = ad.phu.get("ARCHNAME")
            if archname not in arm_dict:
                arm_dict[archname] = {"BLUE": [], "RED": []}
            if "BLUE" in ad.tags:
                arm_dict[archname]["BLUE"].append(ad)
            elif "RED" in ad.tags:
                arm_dict[archname]["RED"].append(ad)
            else:
                log.warning(
                    "No BLUE or RED tag found for %s, skipping this file.", ad.filename
                )

        # Sort streams into red and blue lists
        blue_list = []
        red_list = []
        for archname, arms in arm_dict.items():
            if not arms["BLUE"] or not arms["RED"]:
                msg = f"No BLUE or RED tagged files found for ARCHNAME {archname}"
                log.debug(msg)
                raise ValueError(msg)

            blue_list.extend(arms["BLUE"])
            red_list.extend(arms["RED"])

        log.debug(gt.log_message("primitive", self.myself(), "complete"))
        self.streams["RED"] = red_list
        return blue_list

    def bundleArmStreams(self, adinputs=None, **params):
        """
        Bundle Blue and Red arm AstroData objects.

        This primitive takes the Blue and Red arm streams and combines them
        into multi-extension bundle AstroData objects, reversing the operation
        performed by splitBundle(). Each bundle contains both Blue and Red arms
        as separate extensions.

        Parameters
        ----------
        adinputs : list of AstroData
            List of Blue arm AstroData objects to be combined with
            previously stored Red arm stream in self.streams['RED'].
        suffix : str, optional
            Suffix to append to output filenames. Default is '_reduced'.

        Returns
        -------
        list of AstroData
            List of bundle AstroData objects, each containing Blue and Red
            arm extensions with restored ARCHNAME filenames.

        Notes
        -----
        This primitive requires that the separateArmStreams primitive has
        been run beforehand to populate self.streams['RED'] with the Red
        arm AstroData objects. Each Blue/Red pair must have matching
        ARCHNAME headers to be properly bundled together.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        # Get the red arm stream that was set by separateArmStreams
        if "RED" not in self.streams:
            msg = "RED stream not found. Run separateArmStreams first."
            log.error(msg)
            raise ValueError(msg)

        blue_list = adinputs
        red_list = self.streams["RED"]

        # Create dictionary indexed by ARCHNAME
        blue_dict = {}
        red_dict = {}

        for ad in blue_list:
            archname = ad.phu.get("ARCHNAME")
            if archname is None:
                log.warning("No ARCHNAME found for %s, skipping", ad.filename)
                continue
            blue_dict[archname] = ad

        for ad in red_list:
            archname = ad.phu.get("ARCHNAME")
            if archname is None:
                log.warning("No ARCHNAME found for %s, skipping", ad.filename)
                continue
            red_dict[archname] = ad

        # Check that we have matching pairs
        blue_archnames = set(blue_dict.keys())
        red_archnames = set(red_dict.keys())

        if blue_archnames != red_archnames:
            missing_in_blue = red_archnames - blue_archnames
            missing_in_red = blue_archnames - red_archnames
            if missing_in_blue:
                log.warning("ARCHNAMES in RED but not BLUE: %s", missing_in_blue)
            if missing_in_red:
                log.warning("ARCHNAMES in BLUE but not RED: %s", missing_in_red)

        # Bundle matching pairs
        adoutputs = []
        common_archnames = blue_archnames & red_archnames

        for archname in sorted(common_archnames):
            blue_ad = blue_dict[archname]
            red_ad = red_dict[archname]

            # Create bundle with ARCHNAME as filename
            bundle_ad = deepcopy(blue_ad)
            bundle_ad.filename = archname

            # Update PHU ORIGNAME to archive name
            bundle_ad.phu["ORIGNAME"] = archname

            # Remove ARCHNAME card if it exists (it's now the filename)
            if "ARCHNAME" in bundle_ad.phu:
                del bundle_ad.phu["ARCHNAME"]

            # Append red arm as second extension with all its data
            # This preserves variance, mask, tables, and all other extensions
            bundle_ad.append(red_ad[0])

            # Update name and append to output
            bundle_ad.update_filename(suffix=params["suffix"], strip=True)
            adoutputs.append(bundle_ad)

        log.debug(gt.log_message("primitive", self.myself(), "complete"))
        return adoutputs

    # ========================================================================
    # Private methods
    # ========================================================================

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
        exposuremeter = ad.EXPOSUREMETER.to_pandas(index="Timestamp")
        exposuremeter.index = pd.to_datetime(exposuremeter.index)

        # Reference zeropints
        ref_zp_pc = ad.EXPOSUREMETER.meta["header"]["ZP_PC"]
        ref_zp_frd = ad.EXPOSUREMETER.meta["header"]["ZP_FRD"]

        # 5 min window around exposure for auto-zp determination
        dt1 = TimeDelta(exptime, format="sec")
        dt3 = TimeDelta(300, format="sec")
        utc_end = utc_start + dt1

        # Auto-determine zeropoints if not provided
        n_cutoff = 20
        start_cut = (utc_start - dt3).iso
        end_cut = (utc_end + dt3).iso

        # Fetch context window data unconditionally for both zp determination and plotting
        context_pc_series = exposuremeter.loc[start_cut:end_cut]["Flux PC Channel"]
        context_frd_series = exposuremeter.loc[start_cut:end_cut]["Flux FRD Channel"]
        context_times_pc = Time(context_pc_series.index.values, format="datetime64", scale="utc")
        context_times_frd = Time(context_frd_series.index.values, format="datetime64", scale="utc")
        context_readings_pc = context_pc_series.values.flatten()
        context_readings_frd = context_frd_series.values.flatten()

        if zp_frd == 0.0:
            zp_frd = np.nanmedian(np.sort(context_frd_series)[:n_cutoff])
            log.fullinfo(f"Automatic zeropoint determination for FRD channel: {zp_frd}")
        if zp_pc == 0.0:
            zp_pc = np.nanmedian(np.sort(context_pc_series)[:n_cutoff])
            log.fullinfo(f"Automatic zeropoint determination for PC channel: {zp_pc}")

        # Check if zp are valid
        if np.isnan(zp_frd):
            log.warning("Automatic zeropoint determination for FRD channel failed.")
            zp_frd = 0.0
        if np.isnan(zp_pc):
            log.warning("Automatic zeropoint determination for PC channel failed.")
            zp_pc = 0.0
        if abs(zp_frd - ref_zp_frd) > 0.2 * ref_zp_frd:
            log.warning(
                "Automatic zeropoint determination for FRD, %s, is 20%% off of reference value of %s.",
                zp_frd, ref_zp_frd
            )
        if abs(zp_pc - ref_zp_pc) > 0.2 * ref_zp_pc:
            log.warning(
                "Automatic zeropoint determination for FRD, %s, is 20%% off of reference value of %s.",
                zp_pc, ref_zp_pc
            )

        # Extract readings during exposure
        result_pc = exposuremeter.loc[utc_start.iso : utc_end.iso]["Flux PC Channel"]
        result_frd = exposuremeter.loc[utc_start.iso : utc_end.iso]["Flux FRD Channel"]
        number_pc = result_pc.count()
        number_frd = result_frd.count()

        # Apply zp correction and outlier filtering
        times_pc = Time(result_pc.index.values, format="datetime64", scale="utc")
        readings_pc = result_pc.values.flatten() - zp_pc
        median_pc = medfilt(readings_pc, 3)
        outlier = np.where(np.abs(readings_pc - median_pc) / median_pc > 2)
        readings_pc[outlier] = median_pc[outlier]
        if np.sum(outlier) > 0:
            log.warning("Replaced %s outlier value(s) in PC dataset", len(outlier))
        readings_pc[readings_pc < 0] = 0.0

        times_frd = Time(result_frd.index.values, format="datetime64", scale="utc")
        readings_frd = result_frd.values.flatten() - zp_frd
        median_frd = medfilt(readings_frd, 3)
        outlier = np.where(np.abs(readings_frd - median_pc) / median_frd > 2)
        readings_frd[outlier] = median_frd[outlier]
        if np.sum(outlier) > 0:
            log.warning("Replaced %s outlier value(s) in FRD dataset", len(outlier))
        readings_frd[readings_frd < 0] = 0.0

        # Check for large time gaps in data
        if number_frd > 2 and number_pc > 2:
            gaps = [
                (times_frd[0] - utc_start).sec,
                (utc_end - times_frd[-1]).sec,
                result_frd.index.to_series().diff().dt.total_seconds().fillna(0).max(),
            ]
            maxgap = max(gaps)
            if maxgap > 30:
                log.warning(
                    "%.1f sec gap found in exposuremeter data. Photometric calculations abandoned.",
                    maxgap
                )
            elif maxgap > 10:
                log.warning("%.1f sec gap found in exposuremeter data.", maxgap)
        else:
            times_pc = None
            times_frd = None
            readings_pc = np.nan
            readings_frd = np.nan

        # Calculate statistics
        emeter_stats = {
            "pc": {
                "times": times_pc,
                "readings": readings_pc,
                "stats": {
                    "min": np.nanmin(readings_pc),
                    "max": np.nanmax(readings_pc),
                    "median": np.nanmedian(readings_pc),
                    "std": np.nanstd(readings_pc),
                    "zeropoint": zp_pc,
                },
                "context": {
                    "times": context_times_pc,
                    "readings": context_readings_pc,
                },
            },
            "frd": {
                "times": times_frd,
                "readings": readings_frd,
                "stats": {
                    "min": np.nanmin(readings_frd),
                    "max": np.nanmax(readings_frd),
                    "median": np.nanmedian(readings_frd),
                    "std": np.nanstd(readings_frd),
                    "zeropoint": zp_frd,
                },
                "context": {
                    "times": context_times_frd,
                    "readings": context_readings_frd,
                },
            },
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
        2.0
        * (pars["l"] - spl(1 / mm) * pars["l"])
        * np.cos(pars["theta"])
        * pars["n"]
        / mm
    ) * 1e6


def _fc2min(p, m, etalonwl):
    # residuals are in 'nm' not m/s. Good? bad? Should we normalize?
    return _peak_to_wavelength_spline(m, p) - etalonwl


def _plot_debug(fiber, etalon_peak_data, parameters, residuals):
    """
    Create debug plots for wavelength solution diagnostics.

    Parameters
    ----------
    fiber : int
        Fiber number
    etalon_peak_data : DataFrame
        Peak data with columns: M, WAVELENGTH_BY_THAR, ORDER, CENTER, M_FRACTION
    parameters : lmfit.Parameters
        Etalon parameters used in the fit
    residuals : array
        Residuals in m/s

    Returns
    -------
    fig : matplotlib.figure.Figure
        Debug figure with 6 subplots
    """
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    fig.suptitle(f'Debug Plot - Fiber {fiber}', fontsize=14)

    M = etalon_peak_data["M"].values
    wl = etalon_peak_data["WAVELENGTH_BY_THAR"].values
    orders = etalon_peak_data["ORDER"].values
    center = etalon_peak_data["CENTER"].values

    # Compute predicted wavelength from parameters
    wl_predicted = _peak_to_wavelength_spline(M, parameters)

    # 1. M vs index (color by order)
    ax = axes[0, 0]
    sc = ax.scatter(np.arange(len(M)), M, c=orders, cmap='nipy_spectral', s=2)
    ax.set_xlabel('Index')
    ax.set_ylabel('M (interference order)')
    ax.set_title('M values')
    plt.colorbar(sc, ax=ax, label='Order')

    # 2. WAVELENGTH_BY_THAR vs index (color by order)
    ax = axes[0, 1]
    sc = ax.scatter(np.arange(len(wl)), wl, c=orders, cmap='nipy_spectral', s=2)
    ax.set_xlabel('Index')
    ax.set_ylabel('Wavelength [nm]')
    ax.set_title('WAVELENGTH_BY_THAR')
    plt.colorbar(sc, ax=ax, label='Order')

    # 3. CENTER vs index (color by order)
    ax = axes[1, 0]
    sc = ax.scatter(np.arange(len(center)), center, c=orders, cmap='nipy_spectral', s=2)
    ax.set_xlabel('Index')
    ax.set_ylabel('CENTER [pixels]')
    ax.set_title('CENTER values')
    plt.colorbar(sc, ax=ax, label='Order')

    # 4. Residuals vs wavelength (color by order)
    ax = axes[1, 1]
    sc = ax.scatter(wl, residuals, c=orders, cmap='nipy_spectral', s=2)
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Residuals [m/s]')
    ax.set_title('Residuals vs Wavelength')
    ax.axhline(0, color='k', linestyle='--', alpha=0.5)
    plt.colorbar(sc, ax=ax, label='Order')

    # 5. Predicted vs observed wavelength
    ax = axes[2, 0]
    sc = ax.scatter(wl, wl_predicted, c=orders, cmap='nipy_spectral', s=2)
    ax.set_xlabel('WAVELENGTH_BY_THAR [nm]')
    ax.set_ylabel('Predicted wavelength [nm]')
    ax.set_title('Predicted vs Observed')
    # Add 1:1 line
    lims = [min(wl.min(), wl_predicted.min()), max(wl.max(), wl_predicted.max())]
    ax.plot(lims, lims, 'k--', alpha=0.5)
    plt.colorbar(sc, ax=ax, label='Order')

    # 6. Histogram of residuals
    ax = axes[2, 1]
    valid_res = residuals[~np.isnan(residuals)]
    ax.hist(valid_res, bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Residuals [m/s]')
    ax.set_ylabel('Count')
    ax.set_title(f'Residuals histogram (median={np.nanmedian(residuals):.1f}, std={np.nanstd(residuals):.1f})')
    ax.axvline(np.nanmedian(residuals), color='r', linestyle='--', label='median')
    ax.legend()

    plt.tight_layout()
    return fig