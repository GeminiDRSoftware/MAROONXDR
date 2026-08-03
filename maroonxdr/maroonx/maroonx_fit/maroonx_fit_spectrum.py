"""
Functions describing a single spectrum or a single peak.

The peak model, spectrum evaluation, residual and Jacobian functions
used by the least-squares fits, and the peak finder operating on the
extracted 1D spectrum.
"""
import os
import numpy as np

import scipy.optimize as optimize

from scipy import signal
from scipy.special import erf

from . import get_logger


PLOT_KWARGS = dict(dpi=300, bbox_inches="tight", pad_inches=0.25)

def change_ext(filename, new_ext):
    """
    Change the extension of a filename.

    Parameters
    ----------
    filename : str
        Filename to change.

    new_ext : str
        New extension, without the leading dot.

    Returns
    -------
    str
        Filename with the new extension.
    """
    return os.path.splitext(filename)[0] + "." + new_ext

### POLYMORPHIC DEFINITION FOR SPECTRUM_VALS ###

# Spectrum_vals is polymorphic so that it can be used to
# give the y-values for a given spectrum either by just providing
# a fit object or by providing the parameters and meta parameters.
# Doing this allows us to keep the code cleaner in most places
# as 90% of the time we can just use the MaroonXFit object.x

def peak_val(x, amplitude, center, sigma1, half_width, sigma2):
    """
    Return the values of a peak at the given x values.

    An etalon peak is described by a rectangle step function with a certain
    amplitude, width convolved with a gaussian (with certain sigma). Because
    convolution is numerically intensive, we approximate the convolution with
    a function.

    Parameters
    ----------
    x : ndarray
        x data (pixel).

    amplitude : float
        Peak height.

    center : float
        Peak center position (pixel).

    sigma1 : float
        Sigma of Erf on the left side of the step function (pixel).

    half_width : float
        Half width of the step function (pixel).

    sigma2 : float
        Sigma of Erf on the right side of the step function (pixel).

    Returns
    -------
    ndarray
        Peak values at x.
    """
    arg1 = (x - (center - half_width)) / sigma1
    arg2 = ((center + half_width) - x) / sigma2
    out = erf(arg1) + erf(arg2)
    return 0.5 * amplitude * out

def spectrum_val(x = None, parameters = None, meta_parameters = None, param_obj = None):
    """
    Return the y values of the spectrum model at the given x values.

    Either ``param_obj`` or both ``parameters`` and ``meta_parameters``
    must be given. Each peak is evaluated only inside a window extending
    four sigma beyond the peak edges; outside the windows the model
    equals the offset.

    Parameters
    ----------
    x : ndarray
        x data (pixel).

    parameters : ndarray or None
        Flat fit parameter vector; used together with ``meta_parameters``.

    meta_parameters : MetaParameter or None
        Fit meta parameters; used together with ``parameters``.

    param_obj : Parameter or None
        Parameter object containing the parameters and meta parameters
        of the fit. Takes precedence if given.

    Returns
    -------
    ndarray
        Spectrum values at x.
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
    Return the Jacobian of the spectrum model for the non-center parameters.

    Each peak contributes only inside a window extending six sigma beyond
    the peak edges. The parameters and meta parameters are taken from
    ``fitobj.param_obj``; the ``meta_parameters`` argument is ignored.

    Parameters
    ----------
    bins : ndarray
        x data (pixel) to evaluate the Jacobian at.

    fitobj : MaroonXFit
        Fit object providing the parameters and meta parameters.

    meta_parameters : MetaParameter or None
        Meta parameters. Inherited from legacy but ignored in favour of
        ``fitobj.param_obj.meta_parameters``.

    Returns
    -------
    ndarray
        Jacobian array, one row per x value and one column per
        non-center parameter.
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
    """
    Return the residuals for the polynomial fit.

    Called by ``scipy.optimize.least_squares`` inside
    ``MaroonXFit.fit_polynomials``. The peak-center parameters are taken
    from ``fit_obj.param_obj`` and held fixed. Residuals are clipped to
    three standard deviations around their mean so that outliers in the
    data do not force the solution into a wrong direction; NaNs are
    replaced with zero.

    Parameters
    ----------
    p : ndarray
        Non-center fit parameters proposed by the optimizer.

    fit_obj : MaroonXFit
        Fit object providing the data, fit range, and fixed peak-center
        parameters.

    Returns
    -------
    ndarray
        Clipped residuals of the model against ``fit_obj.data``.
    """

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
    res = np.clip(res,np.nanmean(res)-3*np.nanstd(res),np.nanmean(res)+3*np.nanstd(res))
    return np.nan_to_num(res)

def fit_polynomials_jac(p, fitobj):
    """
    Return the Jacobian for the polynomial fit.

    Called by ``scipy.optimize.least_squares`` inside
    ``MaroonXFit.fit_polynomials``. Updates ``fitobj.param_obj`` in
    place with the proposed parameters before delegating to
    ``spectrum_partial_jacobian``.

    Parameters
    ----------
    p : ndarray
        Non-center fit parameters proposed by the optimizer.

    fitobj : MaroonXFit
        Fit object providing the fit range and parameters.

    Returns
    -------
    ndarray
        Jacobian array from ``spectrum_partial_jacobian``.
    """
    # Get the fixed center parameters
    meta = fitobj.param_obj.meta_parameters
    idx = meta.number_of_peaks
    center_parameters = fitobj.param_obj.parameters[-2 * idx:]
    
    # Reconstruct new parameter vector
    new_parameters = np.concatenate([p, center_parameters])
    fitobj.param_obj.update_parameters(parameters=new_parameters)

    return spectrum_partial_jacobian(fitobj.fitrange, fitobj)

def residual_centers(parameters, x, y, poly_parameters, meta):
    """
    Return the residuals for the peak-center fit.

    Called by ``scipy.optimize.least_squares`` inside
    ``fit_peak_centers``. NaNs are replaced with zero.

    Parameters
    ----------
    parameters : ndarray
        Center and amplitude parameters, two per peak.

    x : ndarray
        x data (pixel).

    y : ndarray
        y data (counts).

    poly_parameters : ndarray
        Polynomial parameters, held fixed.

    meta : MetaParameter
        Fit meta parameters.

    Returns
    -------
    ndarray
        Residuals of the model against ``y``.
    """
    assert len(parameters) == 2 * meta.number_of_peaks, str(p) + str(meta)
    res = spectrum_val(x = x, parameters=np.concatenate([poly_parameters, parameters]), meta_parameters = meta) - y

    return np.nan_to_num(res)


class PeakError(Exception):
    """
    Exception raised when peak finding fails.

    Raised by ``find_peaks``.
    """
    def __str__(self):
        return self.args[0]


def fit_peak_centers(fitrange, data, param_obj, parameter_bounds, iteration = None, fiber = ''):
    """
    Fit the center of each peak individually.

    Each peak's center and amplitude are fitted with
    ``scipy.optimize.least_squares`` inside its bounds window while the
    polynomial parameters are held fixed. The fitted centers and
    amplitudes are stored in place in ``param_obj``. If an individual
    peak fit fails, an error is logged.

    Parameters
    ----------
    fitrange : ndarray
        Domain of the data being fitted (pixel).

    data : ndarray
        1D extracted etalon spectrum.

    param_obj : Parameter
        Object containing the parameters of the fit, along with the meta
        parameters.

    parameter_bounds : ndarray
        Lower and upper parameter bounds of the fit.

    iteration : int or None
        Iteration number, used only in log messages. If None, the log
        messages refer to the guess spectrum instead.

    fiber : str
        Fiber name, used only in log messages.

    Returns
    -------
    list of scipy.optimize.OptimizeResult
        Per-peak least-squares fit results.
    """

    # log = logutils.get_logger(__name__)
    log = get_logger()
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
    Find the local minima and maxima of an etalon spectrum.

    Positive outliers (e.g. cosmic ray hits), values more than 50% above
    the median of the 200 highest samples, are set to NaN in ``data`` in
    place. A Savitzky-Golay filter smooths the data to get rid of noise
    before the extrema search, and the extrema lists are cleaned so that
    minima and maxima alternate. The outermost peaks are not returned to
    avoid problems when fitting the data.

    Parameters
    ----------
    data : ndarray
        1D extracted etalon spectrum. Modified in place: positive
        outliers are set to NaN.

    order : int
        Passed to ``scipy.signal.argrelextrema``. Default is 2.

    savgol_window_length : int
        Window length passed to ``savgol_filter``. Default is 3.

    savgol_polyorder : int
        Polynomial order passed to ``savgol_filter``. Default is 1.

    Returns
    -------
    maxima, minima : tuple of ndarray
        Indices of the local extrema.

    Raises
    ------
    PeakError
        When the cleaned minima and maxima are inconsistent (counts do
        not differ by exactly one or they do not alternate). A debug
        plot is shown before the exception is raised.
    """

    # log = logutils.get_logger(__name__)
    log = get_logger()
    # get rid of noise
    # TODO: maybe improve noise rejection

    # Remove positive outliers (e.g. cosmic ray hits).
    cleanmax = np.sort(data)
    cleanmax = cleanmax[~np.isnan(cleanmax)]
    cleanmax = np.median(cleanmax[-200:])
    nancount = np.count_nonzero(np.isnan(data))
    data[np.array([e>(1.5*cleanmax) if ~np.isnan(e) else False for e in data], dtype=bool) ] = np.nan

    if np.count_nonzero(np.isnan(data)) - nancount > 0:
        log.fullinfo("Removed %d positive outliers 50%% higher than %.3f",\
             np.count_nonzero(np.isnan(data)) - nancount, cleanmax)

    # Clean nan values since they get otherwise smeared out by the filter
    data_clean = data.copy()
    mask = np.isnan(data)
    data_clean[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])

    log.fullinfo(f"Found {np.count_nonzero(np.isnan(data_clean))} nan values in the data")

    # Apply a savgol filter to smooth the data and get rid of noise
    d = signal.savgol_filter(
        data_clean,
        polyorder=savgol_polyorder,
        window_length=savgol_window_length,
        mode="interp",
    )

    log.fullinfo("Applied savgol filter to data")

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

    # -------------------------------------------------
    # Check if there are minima and maxima at the same pixel and remove them
    # minima = np.setdiff1d(minima, maxima)
    # maxima = np.setdiff1d(maxima, minima) # martin

    # Legacy: element-wise comparison (numpy.where(maxima == minima)) removes
    # maxima that coincide with minima at the same index position.
    if len(maxima) == len(minima):
        wrong_maxima = np.where(maxima == minima)
        if np.any(wrong_maxima[0] > 0):
            maxima = np.delete(maxima, wrong_maxima)
    # -------------------------------------------------

    len_minima = len(minima)
    len_maxima = len(maxima)
    # Check that minima and maxima occur alternated, if not, make plots to show error
    if len_minima != len_maxima + 1 or np.any((minima[:-1] > maxima) | (maxima > minima[1:])):
        d_nonnan = d.copy()
        d_nonnan[np.isnan(d)] = 0

        import matplotlib.pyplot as plt
        plt.figure()
        plt.title(f'#Minima: {len(minima)}  #Maxima: {len(maxima)}')
        plt.plot(d)
        plt.plot(maxima, d_nonnan[maxima], 'g+')
        plt.plot(minima, d_nonnan[minima], 'r+')
        plt.hlines(0.25 * np.max(d), 0, len(d))
        plt.show()
        raise PeakError(
            f"Found inconsistent number of peaks: {len(minima)} minima and {len(maxima)} maxima"
            ,(maxima, minima))

    log.fullinfo(f"Found {len_minima} minima and {len_maxima} maxima")
    return maxima, minima

