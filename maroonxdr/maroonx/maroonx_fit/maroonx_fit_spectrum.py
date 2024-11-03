'''
This file contains functions used to describe a single spectrum or a single peak.
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

from scipy import signal
from gempy.utils import logutils
from scipy.special import erf

PLOT_KWARGS = dict(dpi=300, bbox_inches="tight", pad_inches=0.25)

def change_ext(filename, new_ext):
    '''
    Change the extension of a filename
    Parameters
    ----------
    filename: str
        filename to change
    new_ext: str
        new extension
    Returns
    -------
    filename with new extension
    '''
    return os.path.splitext(filename)[0] + "." + new_ext

### POLYMORPHIC DEFINITION FOR SPECTRUM_VALS ###

# Spectrum_vals is polymorphic so that it can be used to
# give the y-values for a given spectrum either by just providing
# a fit object or by providing the parameters and meta parameters.
# Doing this allows us to keep the code cleaner in most places
# as 90% of the time we can just use the MaroonXFit object.x

def peak_val(x, amplitude, center, sigma1, half_width, sigma2):
    """
    Function that returns the values of peaks at a given range for x values.

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
    out = erf(arg1) + erf(arg2)
    return 0.5 * amplitude * out

def spectrum_val(x = None, parameters = None, meta_parameters = None, param_obj = None):
    """
    Function that gives the y values for a given spectrum.

    Parameters
    ----------
    x: np.ndarray
        x data (pixel)
    param_obj: Parameters object containing the parameters and meta parameters
    of the fit

    Returns
    -------
    np.ndarray: spectrum values at x
    """
    if param_obj is not None:
        offset = param_obj.offset
        values = param_obj.eval_polynomials_at_centers()
    else:
        offset = parameters[0]
        centers = parameters[meta_parameters.centers]
        amplitudes = parameters[meta_parameters.amplitudes]
        values = [centers, amplitudes]
        for idx in meta_parameters.polynomials:
            values.append(np.poly1d(parameters[idx])(values[0]))
        values = np.array(values)

    windows = np.zeros(shape=(2, len(values[0])))
    # centers - width -  sigma_left * 4.0
    windows[0] = values[0] - values[4] - 4.0 * values[2]
    # centers + width +  sigma_right * 4.0
    windows[1] = values[0] + values[4] + 4.0 * values[3]

    y = np.full_like(x, offset, dtype=np.float64)
    for values, window in zip(values.T, windows.T):
        l, r = np.searchsorted(x, window)
        p, a, sl, sr, w = values
        y[l:r] += peak_val(x[l:r], a, p, sl, w, sr)
    return y

def spectrum_partial_jacobian(bins, fitobj = None, meta_parameters = None):

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

    param_obj = fitobj.param_obj
    meta_parameters = param_obj.meta_parameters
    values = param_obj.eval_polynomials_at_centers()

    windows = np.zeros(shape=(2, len(values[0])))
    windows[0] = values[0] - values[4] - 6.0 * values[2]
    windows[1] = values[0] + values[4] + 6.0 * values[3]

    rows = meta_parameters.total - 2*meta_parameters.number_of_peaks
    jac = np.zeros(shape=(rows, len(bins)))
    jac[0] = 1.0
    for values, window in zip(values.T, windows.T):
        l, r = np.searchsorted(bins, window)
        x = bins[l:r]
        c, a, sl, sr, w = values
        arg1 = (x - c + w) / sl
        arg2 = (c + w - x) / sr
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

def residual_polynomials(p, fit_obj):
    "Residual function to fit polynomials"

    x = fit_obj.fitrange
    y = fit_obj.data
    meta = fit_obj.param_obj.meta_parameters
    idx = meta.number_of_peaks
    center_parameters = fit_obj.param_obj.parameters[-2 * idx:]
    assert len(center_parameters) == 2*idx, str(p) + str(meta)
    assert len(np.concatenate([p, center_parameters])) == meta.total
    res = spectrum_val(x = x, parameters= np.concatenate([p, center_parameters]),\
        meta_parameters = meta) - y

    # Limit the residuals so that outliers in the data are
    # not forcing the solution into a wrong direction
    p = np.clip(res,np.nanmean(res)-3*np.nanstd(res),np.nanmean(res)+3*np.nanstd(res))
    return np.nan_to_num(p)

def fit_polynomials_jac(_, fitobj):
    "Jacobian for fit_polynomials function"
    return spectrum_partial_jacobian(fitobj.fitrange, fitobj)

def residual_centers(parameters, x, y, poly_parameters, meta):
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
    assert len(parameters) == 2 * meta.number_of_peaks, str(p) + str(meta)
    res = spectrum_val(x = x, parameters=np.concatenate([poly_parameters, parameters]), meta_parameters = meta) - y

    return np.nan_to_num(res)

class PeakError(Exception):
    """
    Exception raised when the fit fails
    """
    def __str__(self):
        return self.args[0]

def fit_peak_centers(fitrange, data, param_obj, parameter_bounds, iteration = None, fiber = ''):
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
    param_obj: Parameters object
        object containing the parameters of the fit, along with the meta parameters
    iteration: int
        Max number of iterations
    fiber: str
        fibre to process, this is only used to generate filenames for plots

    Returns
    -------
    fit_results: list
        list of fit results
    """

    log = logutils.get_logger(__name__)
    # Extract necessary parameters

    parameters = param_obj.parameters
    meta_parameters = param_obj.meta_parameters

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

    # Get the centers and the amplitudes from fit_results so that
    # they can be used to update the parameter obj
    centers = [r.x[0] for r in fit_results]
    amplitudes = [r.x[1] for r in fit_results]

    # Update the parameters object with the new centers and amplitudes.
    # Update done in place, so no need to return
    param_obj.update_parameters(centers=centers, amplitudes=amplitudes)
    if iteration is not None:
        iteration_text = f'{iteration}. vs. {iteration-1}. iteration'
    else:
        iteration_text = 'against guess spectrum'
    log.debug("Mean shift in linecenters for %s (%s): %.2f m/s",\
             fiber, iteration_text, np.nanmedian(np.abs(peak_parameters - centers))*1000.0)
    log.debug ("Bulk shift in linecenters for %s (%s): %.2f m/s",\
         fiber, iteration_text, np.nanmedian(peak_parameters - centers)*1000.0)

    return fit_results

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

    log = logutils.get_logger(__name__)
    # get rid of noise
    # TODO: maybe improve noise rejection

    # Remove positive outliers (e.g. cosmic ray hits).
    cleanmax = np.sort(data)
    cleanmax = cleanmax[~np.isnan(cleanmax)]
    cleanmax = np.median(cleanmax[-200:])
    nancount = np.count_nonzero(np.isnan(data))
    data[np.array([e>(1.5*cleanmax) if ~np.isnan(e) else False for e in data], dtype=bool) ] = np.nan

    if np.count_nonzero(np.isnan(data)) - nancount > 0:
        log.info("Removed %d positive outliers 50%% higher than %.3f",\
             np.count_nonzero(np.isnan(data)) - nancount, cleanmax)

    # Clean nan values since they get otherwise smeared out by the filter
    data_clean = data.copy()
    mask = np.isnan(data)
    data_clean[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])

    log.info(f"Found {np.count_nonzero(np.isnan(data_clean))} nan values in the data")

    # Apply a savgol filter to smooth the data and get rid of noise
    d = signal.savgol_filter(
        data_clean,
        polyorder=savgol_polyorder,
        window_length=savgol_window_length,
        mode="interp",
    )

    log.info("Applied savgol filter to data")

    # use <= and >= for extrema. Sometimes necessary, if signal is 0
    # between peaks over a few pixels. In that case,
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

    # Make sure there is a minima in between each maxima
    k = np.asarray([(np.any(np.logical_and(minima >= maxima[j], minima <= maxima[j + 1])))\
                    for j in range(len(maxima) - 1)])
    missing = np.where(~k)
    for idx in missing:
        minima = np.append(minima,((maxima[idx]+maxima[idx+1])/2).astype('int'))
    minima.sort()

    # Make sure there is a maxima in between each minima
    k = np.asarray([(np.any(np.logical_and(maxima >= minima[j], maxima <= minima[j + 1])))\
                     for j in range(len(minima) - 1)])
    missing = np.where(~k)
    for idx in missing:
        maxima = np.append(maxima,((minima[idx]+minima[idx+1])/2).astype('int'))
    maxima.sort()

    # Check if there are minima and maxima at the same pixel and remove them
    minima = np.setdiff1d(minima, maxima)

    len_minima = len(minima)
    len_maxima = len(maxima)
    # Check that minima and maxima occur alternated, if not, make plots to show error
    if len_minima != len_maxima + 1 or np.any((minima[:-1] > maxima) | (maxima > minima[1:])):
        plt.figure()
        d_nonnan = d.copy()
        d_nonnan[np.isnan(d)] = 0
        plt.title(f'#Minima: {len(minima)}  #Maxima: {len(maxima)}')
        plt.plot(d)
        plt.plot(maxima, d_nonnan[maxima], 'g+')
        plt.plot(minima, d_nonnan[minima], 'r+')
        plt.hlines(0.25 * np.max(d), 0, len(d))
        plt.show()
        raise PeakError(
            f"Found inconsistent number of peaks: {len(minima)} minima and {len(maxima)} maxima"
            ,(maxima, minima))

    log.info(f"Found {len_minima} minima and {len_maxima} maxima")
    return maxima, minima

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
