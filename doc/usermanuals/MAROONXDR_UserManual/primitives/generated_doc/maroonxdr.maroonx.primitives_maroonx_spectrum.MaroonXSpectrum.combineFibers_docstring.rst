
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
adinputs : list of :class:`~astrodata.AstroData`
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
    reaches 10 or clips < max_clips. Default is 5000.

report : bool, optional
    Generate PDF diagnostic report. Default is True.

suffix : str, optional
    Suffix to append to output filenames. Default is empty string.

Returns
-------
list of :class:`~astrodata.AstroData`
    Input frames with the following extensions added (fiber index ``M``
    is 6 by default or 7 when ``symmetric_linefits=True``):

    - ``OPTIMAL_REDUCED_FIBER_M`` : combined flux spectrum with the
      same shape as input fibers.

    - ``OPTIMAL_REDUCED_ERR_M`` : combined error spectrum (inverse of
      summed weights).

    - ``WLS_SIMULTANEOUS_FIBER_6`` or ``WLS_SIMULTANEOUS_SYM_FIBER_7`` :
      wavelength solution copied from fiber 3 (reference fiber).

    - ``REDUCED_ORDERS_FIBER_M`` : list of reduced orders.

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
