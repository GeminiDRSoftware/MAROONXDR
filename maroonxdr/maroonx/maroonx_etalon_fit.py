"""
This file contains helper functions for the fitting of echelle peaks.  These functions by themselves are not primitives,
but are functions that are used in the etalon fitting process.  They are defined here to  prevent primitives_maroonx_echelle.py
from becoming too long and unwieldy.  They are imported into primitives_maroonx_echelle.py and used there.
"""

import logging
import os
import pickle
import time
from collections import namedtuple

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from scipy.special import erf
import scipy.optimize as optimize


log = logging.getLogger(__name__)

PLOT_FORMAT = 'png'
PLOT_DPI = 600

# Add erf to numpy namespace, which allow us to verify our Jacobian
np.erf = erf

def change_ext(filename, new_ext):
    return os.path.splitext(filename)[0] + "." + new_ext

MetaParameterBase = namedtuple(
        "MetaParameterBase",
        [
            "offset",
            "sigma",
            "width",
            "number_of_peaks",
            "total",
            "indices",
            "use_sigma_lr",
        ],   
)

class FitError(Exception):
    """
    Exception raised when the fit fails
    """
    def __str__(self):
        return self.args[0]
    
class MetaParameter(MetaParameterBase):
    """
    namedtuple holding the static meta parameters of the fit and information
    derived from those to  handle the parameter vector

    Note:  Immutable
    """

    #Default values
    OFFSET = 0
    SIGMA_LEFT = 1
    SIGMA_RIGHT = 2
    WIDTH = 3
    CENTERS = 4
    AMPLITUDES = 5

    def __new__(cls, sigma, width, number_of_peaks=0, use_sigma_lr=False):
        offset = 1
        sigma_left = sigma
        sigma_right = sigma if use_sigma_lr else -1
        # Generate a list of indices to split the parameter vector
        split_at = np.cumsum(
            [
                0,
                offset,
                sigma_left + 1,
                sigma_right + 1,
                width + 1,
                number_of_peaks,
                number_of_peaks,
            ]
        )

        # Create a list of slices from the split_at array to actually allow us to split the vector
        indices = [slice(l, r) for l, r in zip(split_at[:-1], split_at[1:])]

        if not use_sigma_lr:  
            # Use the same values for both sigmas- we can do this by making those 2 slices the same
            indices[2] = indices[1]

        total = indices[-1].stop
        self = super().__new__(
            cls, offset, sigma, width, number_of_peaks, total, indices, use_sigma_lr
        )
        return self

    def as_array(self):
        """
        Return as numpy array
        Params
        -------
        self

        Returns
        -------
        A numpy array made from the Meta Parameters
        """
        return np.array(self)

    def change_peaks(self, number_of_peaks):
        """
        Return a copy with changed number of peaks

        Params
        -------
        self
        number_of_peaks: int
        A new number of peaks for the meta parameters

        Returns
        -------
        A copy of the Meta Parameters with the new number of peaks
        """
        return MetaParameter(self.sigma, self.width, number_of_peaks, self.use_sigma_lr)

    def num_polynomial_parameters(self):
        """
        Returns the number of polynomial parameters (incl. offset)
        """
        return self.split_at[-1]

    def __getnewargs__(self):
        """
        Getter function for the argments of the new meta parameters
        """
        return self.sigma, self.width, self.number_of_peaks, self.use_sigma_lr

    @property
    def centers(self):
        """
        Return the indices of the centers of the etalon peaks
        """
        return self.indices[self.CENTERS]

    @property
    def amplitudes(self):
        """
        Return the indices of the amplitudes of the etalon peaks
        """
        return self.indices[self.AMPLITUDES]

    @property
    def polynomials(self):
        """
        Return the indices of the polynomials
        """
        return self.indices[self.SIGMA_LEFT: -2]


def eval_polynomials(x, parameters, meta_parameters):
    """
    Evaluate the polynomials for the given parameters
    Parameters
    ----------
    x: ndarray
        x values to evaluate the polynomials at
    parameters: ndarray
        parameters of the polynomials
    meta_parameters: MetaParameter
        meta parameters of the fit
    Returns
    -------
    values: ndarray
        values of the polynomials at the given x values
    """
    polynomials = meta_parameters.polynomials
    values = np.zeros(shape=(len(polynomials), len(x)))
    for ii, idx in enumerate(polynomials):
        values[ii] = np.poly1d(parameters[idx])(x)
    return values


def get_offset(parameters, meta_parameters):
    """
    Get the offset from the parameters
    """
    return parameters[0] # offset is always the first parameter


def get_centers(parameters, meta_parameters):
    """
    Get the centers of the etalon peaks from the parameters
    """
    return parameters[meta_parameters.centers]


def get_amplitudes(parameters, meta_parameters):
    """
    Get the amplitudes of the etalon peaks from the parameters
    """
    return parameters[meta_parameters.amplitudes]


def eval_polynomials_at_centers(parameters, meta_parameters):
    """
    Evaluate the polynomials at the centers of the fitted etalon peaks
    Parameters
    ----------
    parameters: ndarray
        parameters of the polynomials
    meta_parameters: MetaParameter
        meta parameters of the fit
    Returns
    -------
    ndarray
        values of the polynomials at the centers of the
        fitted etalon peaks
    """
    values = [
        parameters[meta_parameters.centers],
        parameters[meta_parameters.amplitudes],
    ]
    for idx in meta_parameters.polynomials:
        values.append(np.poly1d(parameters[idx])(values[0]))
    return np.array(values)


def concat_parameters(offset, p_sigma_left, p_sigma_right, p_width, peaks, amplitudes, meta_parameters):
    """
    Concatenate a parameter vector into a single numpy array to pass it to
    scipy.least_squares

    Parameters
    ----------
    offset: int
        offset
    p_sigma_left: list of ints
        sigmas to use for the left side
    p_sigma_right: list of ints
        sigmas to use for the right side    
    p_width: list of ints
        coefficients of width polynomial
    peaks: list of ints
        center of etalon peaks
    amplitudes: list of ints
        amplitudes of peaks
    meta_parameters: MetaParameters
        Fit meta parameters

    Returns
    -------
    res : ndarray
        concatenated fit parameters
    """
    if not meta_parameters.use_sigma_lr:
        assert np.all(p_sigma_left == p_sigma_right)
        p_sigma_right = np.array([])

    res = np.concatenate([[offset], p_sigma_left, p_sigma_right, p_width, peaks, amplitudes])
    assert len(res) == meta_parameters.total, f"{len(res)} != {meta_parameters.total}"
    return res


def split_parameters(parameters, meta_parameters):
    """
    Split parameter array into subarrays, reverses concat_parameters

    Parameters
    ----------
    parameters: ndarray
        parameter vector
    meta_parameters: MetaParameters object
        fit meta parameters
    Returns
    ---------
    values: tuple
    (offset, sigma coefficients,
               width coefficient, peak center locations, amplitudes)
    """

    values = [parameters[idx] for idx in meta_parameters.indices]
    centers = values[-2]
    amplitudes = values[-1]
    assert len(centers) == meta_parameters.number_of_peaks, f"{centers} != {meta_parameters.number_of_peaks}"
    assert len(amplitudes) == meta_parameters.number_of_peaks, f"{centers} != {meta_parameters.number_of_peaks}"
    values[0] = np.ndarray.item(values[0])
    return values


def update_parameters(
        parameters,
        meta_parameters,
        offset=None,
        p_sigma_l=None,
        p_sigma_r=None,
        p_width=None,
        centers=None,
        amplitudes=None,
):
    """
    Returns a copy of the parameter vector with updated parameters
    Parameters
    ----------
    parameters: ndarray
        parameter vector
    meta_parameters: MetaParameters object
        fit meta parameters
    offset: int
        offset
    p_sigma_l: list of ints
        sigmas to use for the left side
    p_sigma_r: list of ints
        sigmas to use for the right side
    p_width: list of ints
    Returns
    -------
    p : ndarray
        updated parameter vector
    """
    p = parameters.copy()
    values = [offset, p_sigma_l, p_sigma_r, p_width, centers, amplitudes]
    for idx, v in zip(meta_parameters.indices, values):
        if v is not None:
            p[idx] = v
    return p


class FitResult(object):
    def __init__(
            self,
            x,
            parameters,
            meta_parameters,
            fit_duration,
            iterations,
            polynomial_fit_result,
            peak_fit_results,
    ):
        self.x = x
        self.parameters = parameters
        self.meta_parameters = meta_parameters
        self.fit_duration = fit_duration
        self.iterations = iterations
        self.polynomial_fit_result = polynomial_fit_result
        self.peak_fit_results = peak_fit_results


def find_peaks(data, order=2, savgol_window_length=3, savgol_polyorder=1):
    """
    It performs a savgol filter to smooth the data and get rid of noise.
    The function looks for local minima and maxima. The outer most peaks are
    not returned to avoid problems when fitting the data.

    Parameters
    ----------
    savgol_polyorder: int
        polynomial order passed to savgol_filter
    savgol_window_length: int
        window length passed to savgol_filter
    data: ndarray
        1d extracted etalon spectrum
    order: int
        passed to scipy.signal.argrelextrema
    Returns
    ---------
    maxima, minima: tuple of ndarray
        indices of local extrema: (maxima, minima)
    """
    # get rid of noise
    # TODO: maybe improve noise rejection

    # Remove positive outliers (e.g. cosmic ray hits).
    cleanmax = np.sort(data)
    cleanmax = cleanmax[~np.isnan(cleanmax)]
    cleanmax = np.median(cleanmax[-200:])
    nancount = np.count_nonzero(np.isnan(data))
    data[np.array([e>(1.5*cleanmax) if ~np.isnan(e) else False for e in data], dtype=bool) ] = np.nan

    if np.count_nonzero(np.isnan(data)) - nancount > 0:
        log.info(f'Removed {np.count_nonzero(np.isnan(data)) - nancount} positive outlieres 50% higher than {cleanmax:.3f}')

    # Clean nan values since they get otherwise smeared out by the filter
    data_clean = data.copy()
    mask = np.isnan(data)
    data_clean[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])

    log.info(f'Found {np.count_nonzero(np.isnan(data_clean))} nan values in the data')

    # Apply a savgol filter to smooth the data and get rid of noise
    d = signal.savgol_filter(
        data_clean,
        polyorder=savgol_polyorder,
        window_length=savgol_window_length,
        mode="interp",
    )

    log.info("Applied savgol filter to data")

    # use <= and >= for extrema. Sometimes necessary, if signal is 0 between peaks over a few pixels. In that case,
    # multiple minima are detected, that are cleaned afterwards.
    maxima = signal.argrelextrema(d, np.greater_equal, order=order)[0]
    minima = signal.argrelextrema(d, np.less_equal, order=order)[0]

    def ranges(nums):
        nums = sorted(set(nums))
        gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 1 < e]
        edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
        return list(zip(edges, edges))

    # clean possible double minima / maxima
    minima = np.array([np.mean(r, dtype=int) for r in ranges(minima)]).flatten()
    maxima = np.array([np.mean(r, dtype=int) for r in ranges(maxima)]).flatten()

    # drop incomplete peaks on either end of the spectrum
    left, right = minima[0], minima[-1]
    maxima = maxima[(maxima > left) & (maxima < right)]

    # Make sure there is a minima inbetween each maxima
    k = np.asarray([(np.any(np.logical_and(minima >= maxima[j], minima <= maxima[j + 1]))) for j in range(len(maxima) - 1)])
    missing = np.where(~k)
    for idx in missing:
        minima = np.append(minima,((maxima[idx]+maxima[idx+1])/2).astype('int'))
    minima.sort()

    # Make sure there is a maxima inbetween each minima
    k = np.asarray([(np.any(np.logical_and(maxima >= minima[j], maxima <= minima[j + 1]))) for j in range(len(minima) - 1)])
    missing = np.where(~k)
    for idx in missing:
        maxima = np.append(maxima,((minima[idx]+minima[idx+1])/2).astype('int'))
    maxima.sort()

    # Check if there are minima and maxima at the same pixel and remove them
    minima = np.setdiff1d(minima, maxima)

    # Check that minima and maxima occur alternated, if not, make plots to show error
    if len(minima) != len(maxima) + 1 or np.any((minima[:-1] > maxima) | (maxima > minima[1:])):
        plt.figure()
        d_nonnan = d.copy()
        d_nonnan[np.isnan(d)] = 0
        plt.title(f'#Minima: {len(minima)}  #Maxima: {len(maxima)}')
        plt.plot(d)
        plt.plot(maxima, d_nonnan[maxima], 'g+')
        plt.plot(minima, d_nonnan[minima], 'r+')
        plt.hlines(0.25 * np.max(d), 0, len(d))
        plt.show()
        raise FitError(
             f"Found inconsistent number of peaks: {len(minima)} minima and {len(maxima)} maxima"
             ,(maxima, minima))
    
    log.info(f'Found {len(minima)} minima and {len(maxima)} maxima')
    return maxima, minima

def make_fit_parameters(data, degree_sigma, degree_width, fiber="", plot_path="", use_sigma_lr=True):
    """
    Guess initial parameters for a etalon spectrum fit

    Parameters
    ----------
    data: ndarray
        recorded spectrum
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

    # Given the spectrum, find the peaks
    maxima, minima = find_peaks(data)

    meta_parameters = MetaParameter(degree_sigma, degree_width, len(maxima), use_sigma_lr)
    xx = np.array(range(len(data)))

    # Recalculate maxima based on minima, NOTE: nan values can not be set to zero since 0 for sum of weights not allowed
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

    # Initial guesses
    """
    Fit the polynomial using np.polyfit.  
    Our initial guess is 1.27 + 1.7e-4*maxima for the y-values.
    The degree of the polynomial is meta_parameters.width.
    """
    log.info(f'Initial guess for polynomial fit: 1.27 + 1.7e-4*maxima')
    parameters = concat_parameters(
        offset=np.nanmedian(data[minima]),
        p_sigma_left=np.array([0.0] * meta_parameters.sigma + [0.8]),
        p_sigma_right=np.array([0.0] * meta_parameters.sigma + [0.8]),
        #p_width=np.polyfit(minima, np.diff(minima) / 6.0, meta_parameters.width),
        p_width=np.polyfit(maxima, 1.27+1.7e-4*maxima, meta_parameters.width), #TODO: Replace with something that just buffer adds zeros
        peaks=weighted_pos,
        amplitudes=amplitudes_guess,
        meta_parameters=meta_parameters,
    )

    # Bounds
    parameters_min = concat_parameters(
        -np.inf,
        np.ones(meta_parameters.sigma + 1) * -1,
        np.ones(meta_parameters.sigma + 1) * -1,
        np.ones(meta_parameters.width + 1) * -1,
        minima[:-1].astype(np.float64),
        np.zeros_like(amplitudes_guess),
        meta_parameters=meta_parameters,
    )

    parameters_max = concat_parameters(
        np.inf,
        np.ones(meta_parameters.sigma + 1) * 2.0,
        np.ones(meta_parameters.sigma + 1) * 2.0,
        np.ones(meta_parameters.width + 1) * 2.0,
        minima[1:].astype(np.float64),
        np.full_like(amplitudes_guess, amplitudes_guess*2),
        meta_parameters=meta_parameters,
    )

    # Sanity check
    assert len(parameters_min) == len(parameters), f"{len(parameters_min)} != {len(parameters)}"
    assert len(parameters_max) == len(parameters), f"{len(parameters_max)} != {len(parameters)}"

    parameters_bounds = np.array([parameters_min, parameters_max])

    # Adjust fit range
    fitrange = np.arange(minima[0], minima[-1])

    if plot_path:
        f2 = os.path.join(plot_path, "initial_guess_{}.{}".format(fiber, PLOT_FORMAT))
        plot_fit(fitrange, data[fitrange], parameters, meta_parameters, filename=f2)

    return fitrange, data[fitrange], parameters, parameters_bounds, meta_parameters

def peak(x, amplitude, center, sigma1, half_width, sigma2):
    """
    Function that describes an etalon peak.

    An etalon peak is described by a rectangle step function with a certain
    amplitude, width convolved with a gaussian (with certain sigma). Because
    convolution is numerically intensive, we approximate the convolution with
    a function.

    Parameters
    ----------
    x: np.ndarray
        x data (pixel)
    amplitude: double
        peak height
    center: double
        peak center position (pixel)
    sigma1: double
        sigma of Erf on the left side of the step function (pixel)
    half_width: double
        width of the step function (pixel)
    sigma2: double
        sigma of Erf on the right side of the step function (pixel)

    Returns
    -------
    np.ndarray: peak values at x
    """
    arg1 = (x - (center - half_width)) / sigma1
    arg2 = ((center + half_width) - x) / sigma2
    out = np.erf(arg1) + np.erf(arg2)
    return 0.5 * amplitude * out


def spectrum(x, parameters, meta_parameters):
    """
    Function that describes an full etalon spectrum.

    Parameters
    ----------
    x: np.ndarray
        x data (pixel)
    parameters: np.ndarray
        parameters of the fit
    meta_parameters: MetaParameters
        meta parameters of the fit

    Returns
    -------
    np.ndarray: spectrum values at x
    """

    offset = get_offset(parameters, meta_parameters)
    values = eval_polynomials_at_centers(parameters, meta_parameters)

    y = np.full_like(x, offset, dtype=np.float64)
    for p, a, sl, sr, w, in values.T:
        y += peak(x, a, p, sl, w, sr)
    return y


def spectrum_optimized(x, parameters, meta_parameters):
    """
    Function that describes an full etalon spectrum.

    Parameters
    ----------
    x: np.ndarray
        x data (pixel)
    parameters: np.ndarray
        parameters of the fit
    meta_parameters: MetaParameters
        meta parameters of the fit

    Returns
    -------
    np.ndarray: spectrum values at x
    """

    offset = get_offset(parameters, meta_parameters)
    values = eval_polynomials_at_centers(parameters, meta_parameters)

    windows = np.zeros(shape=(2, len(values[0])))
    # centers - width -  sigma_left * 4.0
    windows[0] = values[0] - values[4] - 4.0 * values[2]
    # centers + width +  sigma_right * 4.0
    windows[1] = values[0] + values[4] + 4.0 * values[3]

    y = np.full_like(x, offset, dtype=np.float64)
    for values, window in zip(values.T, windows.T):
        l, r = np.searchsorted(x, window)
        p, a, sl, sr, w = values
        y[l:r] += peak(x[l:r], a, p, sl, w, sr)
    return y


def spectrum_partial_jacobian_optimized(bins, parameters, meta_parameters):
    """
    Function that describes the jacobian of the spectrum, but ONLY for the
    non center parameters.  Optimized version.

    Parameters
    ----------
    bins: list of ints
        bins
    parameters: np.ndarray
        function parameters
    meta_parameters: MetaParameters object
        meta parameters of the fit
    
    Returns
    ---------
    ndarray: jacobian array
    """

    #print("JAC")

    offset = get_offset(parameters, meta_parameters)
    values = eval_polynomials_at_centers(parameters, meta_parameters)

    windows = np.zeros(shape=(2, len(values[0])))
    # centers - width -  sigma_left * 4.0
    windows[0] = values[0] - values[4] - 6.0 * values[2]
    # centers + width +  sigma_right * 4.0
    windows[1] = values[0] + values[4] + 6.0 * values[3]

    rows = meta_parameters.total - 2*meta_parameters.number_of_peaks
    jac = np.zeros(shape=(rows, len(bins)))
    jac[0] = 1.0
    for values, window in zip(values.T, windows.T):
        l, r = np.searchsorted(bins, window)
        X = bins[l:r]
        c, a, sl, sr, w = values
        arg1 = (X - c + w) / sl
        arg2 = (c + w - X) / sr
        gauss1 = np.exp(-(arg1 ** 2))
        gauss2 = np.exp(-(arg2 ** 2))
        # peaks = 0.5 * (np.erf(arg1) + np.erf(arg2))
        pi_sigmas_l = np.sqrt(np.pi) * sl
        pi_sigmas_r = np.sqrt(np.pi) * sr
        sum_gauss = a * (gauss1 / pi_sigmas_l + gauss2 / pi_sigmas_r)
        sum_arg1_gauss = a * -arg1 * gauss1 / pi_sigmas_l
        sum_arg2_gauss = a * -arg2 * gauss2 / pi_sigmas_r

        idx = 1
        # for ii in range(meta_parameters.amplitude, -1, -1):
        #     jac[idx, l:r] += c ** ii * peaks
        #     idx += 1
        if meta_parameters.use_sigma_lr:
            for ii in range(meta_parameters.sigma, -1, -1):
                jac[idx, l:r] += c ** ii * sum_arg1_gauss
                idx += 1
            for ii in range(meta_parameters.sigma, -1, -1):
                jac[idx, l:r] += c ** ii * sum_arg2_gauss
                idx += 1
        else:
            sum_arg_gauss = a * (-arg1 * gauss1 + -arg2 * gauss2) / pi_sigmas_l
            for ii in range(meta_parameters.sigma, -1, -1):
                jac[idx, l:r] += c ** ii * sum_arg_gauss
                idx += 1
        for ii in range(meta_parameters.width, -1, -1):
            jac[idx, l:r] += c ** ii * sum_gauss
            idx += 1
    return jac.T


def spectrum_partial_jacobian(x, parameters, meta_parameters):
    """
    Function that describes the jacobian of the spectrum, but ONLY for the
    non center parameters.

    Parameters  
    ----------
    x: int
        bins
    parameters: np.ndarray
        function parameters
    meta_parameters: MetaParameters object
        meta parameters of the fit
    
    Returns
    ---------
    jac.T : np.ndarray 
        partial of the jacobian of the spectrum 
    """

    values = eval_polynomials_at_centers(parameters, meta_parameters)
    jac = []
    for c, a, sl, sr, w in values.T:
        arg1 = (x - c + w) / sl
        arg2 = (c + w - x) / sr
        gauss1 = np.exp(-(arg1 ** 2))
        gauss2 = np.exp(-(arg2 ** 2))
        peaks = 0.5 * (np.erf(arg1) + np.erf(arg2))
        pi_sigmas_l = np.sqrt(np.pi) * sl
        pi_sigmas_r = np.sqrt(np.pi) * sr
        sum_gauss = a * (gauss1 / pi_sigmas_l + gauss2 / pi_sigmas_r)

        row = [np.zeros_like(x)]
        # for ii in range(meta_parameters.amplitude, -1, -1):
        #     row.append(c ** ii * peaks)
        if meta_parameters.use_sigma_lr:
            for ii in range(meta_parameters.sigma, -1, -1):
                row.append(c ** ii * a * -arg1 * gauss1 / pi_sigmas_l)
            for ii in range(meta_parameters.sigma, -1, -1):
                row.append(c ** ii * a * -arg2 * gauss2 / pi_sigmas_r)
        else:
            sum_arg_gauss = a * (-arg1 * gauss1 + -arg2 * gauss2)
            pi_sigmas = np.sqrt(np.pi) * sl
            for ii in range(meta_parameters.sigma, -1, -1):
                row.append(c ** ii * sum_arg_gauss / pi_sigmas)
        for ii in range(meta_parameters.width, -1, -1):
            row.append(c ** ii * sum_gauss)

        jac.append(row)
    jac = np.array(jac)
    jac = np.sum(jac, axis=0)
    jac[0] = 1.0  # Offset
    return jac.T


def residual_centers(p, x, y, poly_parameters, meta):
    """
    Residual function to fit the centers.

    Parameters
    ----------
    p: np.ndarray
        parameters
    x: np.ndarray
        x data (pixel)
    y: np.ndarray
        y data (counts)
    poly_parameters: np.ndarray
        polynomial parameters
    meta: MetaParameters object
        meta parameters

    Returns
    -------
    res: np.ndarray
        residual array
    """
    assert len(p) == 2 * meta.number_of_peaks, str(p) + str(meta)
    # res = spectrum_optimized(x, np.concatenate([poly_parameters, p]),
    #                    meta) - y
    res = spectrum_optimized(x, np.concatenate([poly_parameters, p]), meta) - y

    return np.nan_to_num(res)


def fit_peak_centers(fitrange, data, parameters, parameter_bounds, meta_parameters, iteration = None, fiber = ''):
    """
    Function to fit the center of each peak individually

    Parameters
    ----------
    fitrange: np.ndarray
        domain of the data being fitted
    data: np.ndarray
        1D extracted etalon spectrum
    parameters: np.ndarray
        parameters of the fit
    parameter_bounds: np.ndarray
        parameter bounds of the fit
    meta_parameters: MetaParameters object
        meta parameters of the fit
    iteration: int  
        Max number of iterations
    fiber: str
        fibre to process, this is only used to generate filenames for plots
    
    Returns
    -------
    result: np.ndarray
        fitted parameters with the peak centers fixed
    fit_results: list
        list of fit results
    """

    # Extract necessary parameters
    idx = meta_parameters.number_of_peaks
    meta_parameters_peak = meta_parameters.change_peaks(1)
    poly_parameters = parameters[: -2 * idx]
    peak_parameters = parameters[-2 * idx: -idx]
    amplitude_parameters = parameters[-idx:]
    peak_parameter_bounds = parameter_bounds.T[-2 * idx: -idx]
    amplitude_parameter_bounds = parameter_bounds.T[-idx:]

    # Fit the peak centers
    fit_results = []
    for p, amplitude, peakbounds, amplitudebounds in zip(
        peak_parameters,
        amplitude_parameters,
        peak_parameter_bounds,
        amplitude_parameter_bounds,
    ):
        l, r = np.searchsorted(fitrange, peakbounds).astype(int)
        x_peak = fitrange[l:r]
        data_peak = data[l:r]
        res = optimize.least_squares(
            residual_centers,
            [p, amplitude],
            args=(x_peak, data_peak, poly_parameters, meta_parameters_peak),
            bounds=(
                [peakbounds[0], amplitudebounds[0]],
                [peakbounds[1], amplitudebounds[1]+0.001],
            ),
            x_scale='jac',
        )
        if not res.success:
            log.error(f"Failed fitting peak at {p}.")
        fit_results.append(res)

    centers = [r.x[0] for r in fit_results]
    amplitudes = [r.x[1] for r in fit_results]
    result = update_parameters(
        parameters, meta_parameters, centers=centers, amplitudes=amplitudes
    )
    if iteration is not None:
        iteration_text = f'{iteration}. vs. {iteration-1}. iteration'
    else:
        iteration_text = 'against guess spectrum'
    log.debug(f'Mean shift in linecenters for {fiber} ({iteration_text}):  {np.nanmedian(np.abs(peak_parameters - centers))*1000.0:.2f} m/s')
    log.debug(f'Bulk shift in linecenters for {fiber} ({iteration_text}):  {np.nanmedian(peak_parameters - centers) * 1000.0:.2f} m/s')
    return result, fit_results


def residual_polynomials(p, x, y, center_parameters, meta):
    
    "Residual function to fit polynomials"
    assert len(center_parameters) == 2*meta.number_of_peaks, str(p) + str(meta)
    assert len(np.concatenate([p, center_parameters])) == meta.total
    res = spectrum_optimized(x, np.concatenate([p, center_parameters]), meta) - y

    # Limit the residuals so that outliers in the data are not forcing the solution into a wrong direction
    res = np.clip(res,np.nanmean(res)-3*np.nanstd(res),np.nanmean(res)+3*np.nanstd(res))

    return np.nan_to_num(res)


def fit_polynomials_jac(p, x, _, center_parameters, meta):
    "Jacobian for fit_polynomials function"
    return spectrum_partial_jacobian_optimized(
        x, np.concatenate([p, center_parameters]), meta
    )


def fit_polynomials(fitrange, data, parameters, parameter_bounds, meta_parameters, **least_square_kw):
    """
    Fit the the polynomial parameters while holding the peak centers.

    Parameters
    ----------
    fitrange: np.ndarray
        pixel range of the data
    data: np.ndarray
        1D extracted etalon spectrum
    parameters: np.ndarray
        parameters of the fit
    parameter_bounds: np.ndarray
        parameter bounds of the fit
    meta_parameters: MetaParameters object
        meta parameters of the fit
    
    Returns
    -------
    result: np.ndarray
        fitted parameters with the peak centers fixed
    """

    kw_args = dict(xtol=1e-15, ftol=1e-12)
    kw_args.update(least_square_kw)

    idx = meta_parameters.number_of_peaks
    result = optimize.least_squares(
        residual_polynomials,
        parameters[:-2*idx],
        args=(fitrange, data, parameters[-2*idx:], meta_parameters),
        bounds=parameter_bounds[:, :-2*idx],
        jac=fit_polynomials_jac,
        x_scale='jac',
        **kw_args
    )

    if not result.success:
        raise FitError("Error fitting polynomials.", result)
    # The values for the polynomials for each peak should be larger than 0
    evaluated_parameters = eval_polynomials_at_centers(parameters, meta_parameters)
    if np.any(evaluated_parameters < 0):
        raise FitError("Error fitting polynomials: Amplitudes, sigma or width was < 0", result)
    return np.concatenate([result.x, parameters[-2*idx:]]), result


def plot_peaks(data, minima, maxima, ax=None, filename=None):
    """
    Plot the result of find_peaks

    Parameters
    ----------
    data: np.array like
        1d extracted etalon spectrum
    minima: np.array like
        indices of local minima
    maxima: np.array like
        indices of local maxima
    ax: matplotlib.Axis
        Axis to use for plotting
    filename: str   
        file to save plot

    Returns
    -------
    Plots added to the filename
    """

    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(data, ".--", color="grey")
    data_nonnan = data.copy()
    data_nonnan[np.isnan(data)] = 0
    ax.plot(minima, data_nonnan[minima], ".", color="blue", label="Minima")
    ax.plot(maxima, data_nonnan[maxima], ".", color="red", label="Maxima")
    ax.set_title("Found peaks")

    if filename is not None:
        plt.savefig(filename, **PLOT_KWARGS)
        np.savez_compressed(
            change_ext(filename, "npz"), data=data, minima=minima, maxima=maxima
        )
        plt.clf()
        plt.cla()


def plot_fit(fitrange, data, parameters, meta_parameters, ax1=None, ax2=None, filename=None):
    """
    Plot the result of fit_polynomials

    Parameters
    ----------
    fitrange: np.array like
        x data (pixel) - domain of the data being fitted
    data: np.array like
        1d extracted etalon spectrum
    maxima: list / np.nd.array
        indices of local maxima
    parameters: np.array
        parameters that describe etalon spectrum
    meta_parameters: MetaParameters object
        meta parameters of the fit
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
    # Extract necessary parameters
    fit = spectrum(fitrange, parameters, meta_parameters)
    centers = get_centers(parameters, meta_parameters)
    polynomials = eval_polynomials(fitrange, parameters, meta_parameters)

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
            change_ext(filename, "npz"),
            data=data,
            fit=fit,
            centers=centers,
            x=fitrange,
            polynomials=polynomials,
        )
        plt.clf()
        plt.cla()


def iterative_fit(
        input_spectrum,
        degree_sigma,
        degree_width,
        iterations=3,
        guess_spectrum=None,
        fiber="",
        plot_path="",
        use_sigma_lr=True,
        show_plots=False
):
    """Fits model to etalon spectrum.

    QUICK-FIX (2021-09-08, AS): Replace (not supported) guess parameters from guess file with guess spectrum
    from this file to run guess on the data. Useful not as a 'speed-up' but for potentially bad/contaminated
    data that need a good first guess

    TODO: possible refinements:
        - Evaluate the cost function or the change in peak center positions
          as stop criterion.
        - An initial_guess could be used. Typically, we will have a set of data
          that is very similar to each other (e.g. etalon spectra of one day).
          The initial_guess can be used to speed up fitting as long as the
          changes in between spectra are small.

    TODO: at the moment there is only one global function for sigma and sigma1
    and sigma2 is forced to be identical. However, it is more precise to
    describe left and right wings of the peaks with seperate sigmas. -> make
    two global functions to make seperate functional describtion of sigma for
    left and right

    Parameters
    ----------
    spectrum_: np.ndarray
        recorded etalon spectrum
    degree_sigma: int
        degree for the polynomial of the peak sigmas
    degree_width: int
        degree for the polynomial of the peak widths
    iterations: int
        Number of iterations to fit
    guess_spectrum: np.ndarray
        Extracted etalon spectrum from a previous fit (see TODO)
    fibre: str
        fibre to process, this is only used to generate filenames for
        debug plots and for log messages
    plot_path: str
        Generate debug plots and save them at the given location
    
    Returns
    -------
    FitResult: namedtuple
    """

    t0 = time.time()
    # If there is a guess spectrum, use this for the first iteration of fitting data
    if guess_spectrum is not None:
        fit_spectrum = guess_spectrum
        log.info("Use guess spectrum for initial fit")
    else:
        # Otherwise use the input spectrum for the first iteration
        fit_spectrum = input_spectrum

    # Create fit parameters for initial fit
    fitrange, data, parameters, bounds, meta_parameters = make_fit_parameters(
        data=fit_spectrum,
        degree_sigma=degree_sigma,
        degree_width=degree_width,
        fiber=fiber,
        plot_path=plot_path,
        use_sigma_lr=use_sigma_lr,
    )

    # If provided with the guess spectrum, fit the peak centers using the parameters from the guess file
    if guess_spectrum is not None:
        data = input_spectrum[fitrange]
        parameters, fit_results = fit_peak_centers(
            fitrange, data, parameters, bounds, meta_parameters
        )

    shw = eval_polynomials(fitrange, parameters, meta_parameters)[2, :]
    log.info(f'Slit half-width on {fiber} (0. iteration): {shw[0]:.3f} - {shw[-1]:.3f}')
    amplitudes_sorted = np.sort(eval_polynomials_at_centers(parameters, meta_parameters)[1, :])
    log.info(f'Flat-relative amplitudes on {fiber} (0. iteration): min: {amplitudes_sorted[10]:.2f}, \
                max: {amplitudes_sorted[-10]:.2f}, median: {np.median(amplitudes_sorted):.2f}')

    fit = spectrum(fitrange, parameters, meta_parameters)
    residuals = np.abs((data-fit)/fit)
    bad = np.where(residuals > 2)
    residuals = np.clip(residuals,0,2)
    residuals[np.where(residuals > 5*np.nanstd(residuals))] = np.nan
    resis = np.nansum(residuals)
    samples = np.sum(~np.isnan(residuals))
    log.info(f'Fit residuals on {fiber} (0. iteration): {resis:.1f} ({resis/samples:.3f})')

    log.info(f"Clip {np.count_nonzero(bad)} datapoints as outliers when calculating residuals")
    data[bad] = np.nan

    log.info(f'Degree sigma:{degree_sigma}, degree width: {degree_width}')
    log.debug("Created fit parameters for %s in %is", fiber, time.time() - t0)

    ############################################################
    # CODE BLOCK
    #
    # if initial_parameters is not None:
    #     initial_parameters, initial_meta_parameters = initial_parameters
    #     logger.info("Use initial parameter for %s", fiber)
    #     if initial_meta_parameters == meta_parameters:
    #         logger.error("Initial parameters are not yet implemented")
    #         # TODO if the initial_parameters are used, we need to ensure that
    #         # the number of peaks and the polynomial degrees are compatible.
    #         # Further we # should keep the bounds detected by
    #         # make_fit_parameters, because the peaks might shift
    #         #
    #         # parameters = initial_parameters
    #         # iterations -= 1
    #     else:
    #         logger.error(
    #             "Initial parameters don't match %s %s: ignoring them",
    #             initial_meta_parameters,
    #             meta_parameters,
    #         )
    #     # raise NotImplementedError("Initial parameters are not supported yet")
    ############################################################

    if show_plots:
        plot_fit(fitrange, data, parameters, meta_parameters)
        plt.show()

    t0 = time.time()

    ############################################################
    # CODE BLOCK
    #
    # parameters, _ = fit_peak_centers(x, data, parameters, bounds, meta_parameters)
    # logger.debug(
    #     "Fitted centers (1. iteration) on %s in: %.2fs", fiber, time.time() - t0
    # )
    #
    # if show_plots:
    #     plot_fit(x, data, parameters, meta_parameters)
    #     plt.show()
    # if plot_path:
    #     f1 = os.path.join(plot_path, "initial_{}.{}".format(fiber, PLOT_FORMAT))
    #     plot_fit(x, data, parameters, meta_parameters, filename=f1)
    ############################################################

    p = [parameters.copy()]
    old_resis_norm = 10
    converged = False
    for ii in range(1, iterations):
        t_ii = time.time()
        parameters, resf = fit_polynomials(fitrange, data, parameters, bounds, meta_parameters)
        log.debug(
            "Fitted polynomials (%i. iteration) on %s in: %.2fs",
            ii + 1,
            fiber,
            time.time() - t_ii,
        )
        shw = eval_polynomials(fitrange, parameters, meta_parameters)[2, :]
        log.debug(f'Slit half-width on {fiber} ({ii}. iteration): {shw[0]:.3f} - {shw[-1]:.3f}')
        amplitudes_sorted = np.sort(eval_polynomials_at_centers(parameters,meta_parameters)[1,:])
        log.debug(f'Flat-relative amplitudes on {fiber} ({ii}. iteration): min: {amplitudes_sorted[10]:.2f}, \
                    max: {amplitudes_sorted[-10]:.2f}, median: {np.median(amplitudes_sorted):.2f}')
        p.append(parameters.copy())

        t_ii = time.time()
        parameters, fit_results = fit_peak_centers(
            fitrange, data, parameters, bounds, meta_parameters, iteration=ii, fiber=fiber
        )

        log.debug(
            "Fitted centers (%i. iteration) on %s in: %.2fs",
            ii + 1,
            fiber,
            time.time() - t_ii,
        )

        fit = spectrum(fitrange, parameters, meta_parameters)
        residuals = np.abs((data - fit) / fit)
        residuals = np.clip(residuals, 0, 1)
        residuals[np.where(residuals > 3 * np.nanstd(residuals))] = np.nan
        resis = np.nansum(residuals)
        samples = np.sum(~np.isnan(residuals))
        resis_norm = resis / samples
        log.debug(f'Fit residuals on {fiber} ({ii}. iteration): {resis:.1f} ({resis_norm:.3f})')
        if np.abs(resis_norm - old_resis_norm) < 0.05 * old_resis_norm:
            # If new residuals differ less than 5% from old residuals, stop iterating
            converged = True
            if resis_norm < old_resis_norm:
                # if we still see in improvement, update results, then break
                ii = ii + 1
                p.append(parameters.copy())
            break
        else:
            p.append(parameters.copy())
        old_resis_norm = resis_norm

        if plot_path:
            f1 = os.path.join(plot_path, "fit_results_{}_{}.pkl".format(ii, fiber))
            with open(f1, "bw") as outfile:
                pickle.dump((resf, fit_results), outfile)

    if converged:
        log.info(f'No more improvement for {fiber} after {ii-1}. iteration, finished')
    else:
        log.warning(f'Global fit did not converge for {fiber} after {ii+1} iterations')

    if show_plots:
        plot_fit(fitrange, data, parameters, meta_parameters)
        plt.show()

    if plot_path:
        f1 = os.path.join(plot_path, "fit_{}.{}".format(fiber, PLOT_FORMAT))
        plot_fit(fitrange, data, p[-1], meta_parameters, filename=f1)
        f2 = os.path.join(plot_path, "iterations_{}.{}".format(fiber, "npz"))
        np.savez_compressed(f2, *p, x=fitrange, spectrum=spectrum, data=data, bounds=bounds)
        f3 = os.path.join(plot_path, "meta_parameters_{}.{}".format(fiber, "pkl"))
        with open(f3, "bw") as outfile:
            pickle.dump(meta_parameters, outfile)

    t_total = time.time() - t0
    log.info(f"Fitting {fiber} done in {t_total :.2f} s", fiber, t_total)
    return FitResult(fitrange, p[-1], meta_parameters, t_total, iterations, resf, fit_results)


        

