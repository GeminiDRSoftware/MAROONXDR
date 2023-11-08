'''
This function contains the MaroonXFit class,
which is used to describe the parameters of a fit to a MaroonX spectrum.
It also contains other, smaller, related classes.  Functions in maroonx_fit_spectrum are
overloaded here to use the MaroonXFit object.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from gempy.utils import logutils

from . import maroonx_fit_spectrum as spectrum
from .maroonx_fit_parameters import Parameter, MetaParameter


PLOT_KWARGS = dict(dpi=300, bbox_inches="tight", pad_inches=0.25)

class FitError(Exception):
    '''
    Exception raised when the fit fails
    '''
    def __str__(self):
        return self.args[0]


class FitResult(object):
    '''
    This class is used to describe the result of a fit to a MaroonX spectrum.
    Contains the following properties:
    fit_obj: MaroonXFit object
        MaroonXFit object that was fitted
    fit_duration: float
        time it took to fit the spectrum
    iterations: int
        number of iterations it took to fit the spectrum
    polynomial_fit_result: np.array
        result of the polynomial fit
    peak_fit_results: np.array
        result of the peak fit
    '''
    def __init__(
            self,
            fit_obj,
            fit_duration,
            iterations,
            polynomial_fit_result,
            peak_fit_results):
        self.fit_obj = fit_obj
        self.fit_duration = fit_duration
        self.iterations = iterations
        self.polynomial_fit_result = polynomial_fit_result
        self.peak_fit_results = peak_fit_results

class MaroonXFit(object):
    '''
    This class is used to describe the parameters of a fit to a MaroonX spectrum.
    Contains the following properties:
    fitrange: np.array like
        x data (pixel) - domain of the data being fitted
    data: np.array like
        1d extracted etalon spectrum
    param_obj: Parameter object
        parameters that describe etalon spectrum
    parameters_bounds: np.array
        parameter bounds of the fit
    fiber : str
        fiber to process, this is only used to generate filenames for
        debug plots and for log messages
    plot_path: str
        Generate debug plots and save them at the given location
    '''
    def __init__(self, data, degree_sigma=None, degree_width=None,\
                  fiber="", plot_path="", use_sigma_lr=True):
        """
        Guess initial parameters for a etalon spectrum fit

        Parameters
        ----------
        data: ndarray
            recorded 1D spectrum
        degree_sigma: tuple
            degree of the sigma polynomial
        degree_width: int
            degree of the half width polynomial
        fiber: str
            fibre to process, this is only used to generate filenames for
            debug plots and for log messages
        plot_path: str
            Generate debug plots and save them at the given location
        use_sigma_lr: bool
            Use different polynomial for sigma at left and
            right flank
        Returns
        --------
        tuple: (initial parameters, parameter boundaries, fit meta parameters)
        The tuple contains the initial parameters for the fit, and the boundaries
        """

        self.log = logutils.get_logger(__name__)
        log = self.log
        # Given the spectrum, find the peaks
        maxima, minima = spectrum.find_peaks(data)

        meta_parameters = MetaParameter(degree_sigma, degree_width, len(maxima), use_sigma_lr)
        xx = np.array(range(len(data)))

        # Recalculate maxima based on minima,
        # NOTE: nan values can not be set to zero since 0 for sum of weights not allowed
        weighted_pos = np.array(
            [
                np.average(
                    xx[minima[idx_minima]: minima[idx_minima + 1]],
                    weights=np.nan_to_num(
                        data[minima[idx_minima]: minima[idx_minima + 1]],nan=1e-5),
                )
                for idx_minima in range(len(minima) - 1)
            ]
        )
        log.info("Recalculated maxima based on minima")

        # Find bad amplitudes by checking where we have NaNs at the maxima
        amplitudes_guess = data[maxima]
        bad = np.where(np.isnan(amplitudes_guess))
        if bad is not None:
            # Replace bad amplitudes by taking means without NaNs
            for b in bad:
                amplitudes_guess[b] = np.nanmean(amplitudes_guess)


        # Fit the polynomial using np.polyfit.
        # Our initial guess is 1.27 + 1.7e-4*maxima for the y-values.
        # The degree of the polynomial is meta_parameters.width.


        log.info("Initial guess for polynomial fit: 1.27 + 1.7e-4*maxima")

        param_obj = Parameter(
            offset=np.nanmedian(data[minima]),
            p_sigma_left=np.array([0.0] * meta_parameters.sigma + [0.8]),
            p_sigma_right=np.array([0.0] * meta_parameters.sigma + [0.8]),
            p_width=np.polyfit(maxima, 1.27+1.7e-4*maxima, meta_parameters.width),
            #TODO: Replace with something that just buffer adds zeros
            peaks=weighted_pos,
            amplitudes=amplitudes_guess,
            meta_parameters=meta_parameters)

        # Bounds
        minimum = Parameter(
            -np.inf,
            np.ones(meta_parameters.sigma + 1) * -1,
            np.ones(meta_parameters.sigma + 1) * -1,
            np.ones(meta_parameters.width + 1) * -1,
            minima[:-1].astype(np.float64),
            np.zeros_like(amplitudes_guess),
            meta_parameters=meta_parameters)

        maximum = Parameter(
            np.inf,
            np.ones(meta_parameters.sigma + 1) * 2.0,
            np.ones(meta_parameters.sigma + 1) * 2.0,
            np.ones(meta_parameters.width + 1) * 2.0,
            minima[1:].astype(np.float64),
            np.full_like(amplitudes_guess, amplitudes_guess*2),
            meta_parameters=meta_parameters)

        parameters_min = minimum.parameters
        parameters_max = maximum.parameters
        parameters = param_obj.parameters
        # Sanity check
        assert len(parameters_min) == len(parameters),\
        f"{len(parameters_min)} != {len(parameters)}"
        assert len(parameters_max) == len(parameters),\
            f"{len(parameters_max)} != {len(parameters)}"

        parameters_bounds = np.array([parameters_min, parameters_max])

        # Adjust fit range
        fitrange = np.arange(minima[0], minima[-1])

        self.fitrange = fitrange
        self.data = data[fitrange]
        self.param_obj = param_obj
        self.parameters_bounds = parameters_bounds
        self.fiber = fiber
        self.plot_path = plot_path


    def fit_polynomials(self, **least_square_kw):
        """
        Fit the the polynomial parameters while holding the peak centers.

        Parameters
        ----------
        Maroon_X_Fit object: Maroon_X_Fit object that is being fitted

        Returns
        -------
        result: np.ndarray
            fitted parameters with the peak centers fixed
        """

        kw_args = dict(xtol=1e-15, ftol=1e-12)
        kw_args.update(least_square_kw)

        # Extract necessary parameters
        param_obj = self.param_obj
        parameter_bounds = self.parameters_bounds

        parameters = param_obj.parameters
        meta_parameters = param_obj.meta_parameters

        idx = meta_parameters.number_of_peaks

        #Update the parameters for the fit
        result = optimize.least_squares(
        spectrum.residual_polynomials,
        param_obj.parameters[:-2*idx],
        args = (self,),
        bounds=parameter_bounds[:, :-2*idx],
        jac= spectrum.fit_polynomials_jac,
        x_scale='jac',
        **kw_args)

        if not result.success:
            raise FitError("Error fitting polynomials.", result)
        # The values for the polynomials for each peak should be larger than 0
        evaluated_parameters = param_obj.eval_polynomials_at_centers()
        if np.any(evaluated_parameters < 0):
            raise FitError("Error fitting polynomials: Amplitudes, sigma or width was < 0", result)

        parameters = np.concatenate([result.x, parameters[-2*idx:]])
        updated_param_obj = param_obj.update_parameters(parameters=parameters)
        self.update_maroonxfit(param_obj=updated_param_obj)
        return result

    def plot_fit(self, ax1=None, ax2=None, filename=None):
        """
        Plot the result of fit_polynomials

        Parameters
        ----------
        Maroon_X_Fit object: Maroon_X_Fit object that is being fitted
        ax1: matplotlib.Axis
            Axis to use for plotting peaks (optional)
        ax2: matplotlib.Axis
            Axis to use for plotting ?? (optional)
        filename: str
            file to save plot (optional)

        Returns
        -------
        Plots added to file with filename
        """

        data = self.data
        fitrange = self.fitrange
        # Extract necessary parameters
        fit = self.spectrum_val()
        centers = self.centers()
        polynomials = self.eval_polynomials()

        # Plot
        if ax1 is None:
            fig = plt.figure(constrained_layout=True)
            gs = fig.add_gridspec(3, 1)
            ax1 = fig.add_subplot(gs[0, 0])
            ax3 = fig.add_subplot(gs[1, 0],sharex=ax1)
            ax2 = fig.add_subplot(gs[2, 0],sharex=ax1)

        # Calculate residuals
        residuals = np.abs((data-fit)/fit)
        residuals = np.clip(residuals,0,2)
        residuals[np.where(residuals > 5*np.nanstd(residuals))] = np.nan
        resis = np.nansum(residuals)

        # Create the fit plot
        ax1.set_title("Fit")
        ax1.plot(fitrange, data, label="data")
        ax1.plot(centers, np.interp(centers, fitrange, data), "x", color="grey", label="peaks")
        ax1.plot(fitrange, fit, color="blue", label="fit")
        data_sorted = np.sort(data)
        ax1.set_ylim(-0.01, 1.1*data_sorted[-50])
        ax1.legend()

        # Create the residual plot
        ax3.set_title(f'Residuals ({resis:.1f})')
        ax3.plot(fitrange, data-fit, label="Residuals")
        dmf = data-fit
        dmf[np.where(np.isnan(residuals))] = np.nan
        ax3_yl = np.nanmin(dmf)
        ax3_yu = np.nanmax(dmf)
        ax3.set_ylim(ax3_yl,ax3_yu)

        # Create the plot with the polynomial values
        ax2.set_title("Polynomial values")
        labels = ["amplitude", "sigma_l", "sigma_r", "box half width"]
        labels = ["sigma_l", "sigma_r", "width"]
        for label, values in zip(labels, polynomials):
            ax2.plot(fitrange, values, label=label)
        ax2.legend()

        if filename is not None:
            plt.savefig(filename, **PLOT_KWARGS)
            np.savez_compressed(
                spectrum.change_ext(filename, "npz"),
                data=data,
                fit=fit,
                centers=centers,
                x=fitrange,
                polynomials=polynomials,
            )
            plt.clf()
            plt.cla()

    def update_maroonxfit(self, fitrange = None, data = None, param_obj = None, param_bounds = None):
        """
        Update the MaroonXFit object with new data and parameters

        Parameters
        ----------
        Maroon_X_Fit object: Maroon_X_Fit object that is being fitted
        fitrange: np.array like or None
            x data (pixel) - domain of the data being fitted.  If None, the
            fitrange of the object is not updated
        data: np.array like or None
            1d extracted etalon spectrum. If None, the data of the object is
            not updated
        param_obj: Parameter object or None
            parameters that describe etalon spectrum. If None, the parameters
            of the object is not updated
        param_bounds: np.array or None
            parameter bounds of the fit. If None, the parameter bounds of the
            object is not updated
        Returns
        -------
        None
        """
        if fitrange is not None:
            self.fitrange = fitrange
        if data is not None:
            self.data = data
        if param_obj is not None:
            self.param_obj = param_obj.copy()
        if param_bounds is not None:
            self.parameters_bounds = param_bounds


    # The following functions are overloaded versions of functions in maroonx_fit_spectrum.py

    def centers(self):
        """
        Overloads the centers function in maroonx_fit_parameters.py.
        Function to calculate the centers of the fitted etalon peaks

        Parameters
        ----------
        Maroon_X_Fit object

        Returns
        -------
        centers: ndarray
            centers of the fitted etalon peaks
        """
        param_obj = self.param_obj
        return param_obj.centers

    def spectrum_val(self):
        """
        Overloads the spectrum_val function in maroonx_fit_spectrum.py.
        Function to calculate the spectrum values for the given Maroon_X_Fit object

        Parameters
        ----------
        Maroon_X_Fit object

        Returns
        -------
        fit: ndarray
            values of the fit for x values in fitrange
        """
        x = self.fitrange
        param_obj = self.param_obj
        return spectrum.spectrum_val(x = x, param_obj = param_obj)

    def eval_polynomials(self):
        """
        Overloads the eval_polynomials function in maroonx_fit_parameters.py
        Evaluates the polynomials for the given Maroon_X_Fit object using the
        parameters and fitrange provided in the object

        Parameters
        ----------
        Maroon_X_Fit object

        Returns
        -------
        values: ndarray
            values of the polynomials for x values in fitrange
        """
        param_obj = self.param_obj
        fitrange = self.fitrange
        return param_obj.eval_polynomials(fitrange)

    def eval_polynomials_at_centers(self):
        """
        Overloads the eval_polynomials_at_centers function in maroonx_fit_parameters.py.
        Evaluate the polynomials at the centers of the fitted etalon peaks.
        Parameters
        ----------
        Maroon_X_Fit object

        Returns
        -------
        ndarray
            values of the polynomials at the centers of the
            fitted etalon peaks
        """
        param_obj = self.param_obj
        return param_obj.eval_polynomials_at_centers()

    def fit_peak_centers(self, iteration = None, fiber = None):
        """
        Overloads the fit_peak_centers function in maroonx_fit_spectrum.py.
        Function to fit the center of each peak individually
        Parameters
        ----------
        Maroon_X_Fit object

        Returns
        -------
        result: np.ndarray
            fitted parameters with the peak centers fixed
        """
        fitrange = self.fitrange
        data = self.data
        param_obj = self.param_obj
        parameters_bounds = self.parameters_bounds

        fit_results = spectrum.fit_peak_centers(fitrange, data, param_obj,\
                     parameters_bounds, iteration, fiber)
        return fit_results
