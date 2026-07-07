
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
adinputs : list of :class:`~astrodata.AstroData`
    Input AstroData objects containing 1D extracted science spectra with
    PEAKS and POLY extensions from getPeaksAndPolynomials.

wavecal : str or :class:`~astrodata.AstroData`, optional
    Corresponding etalon calibration file with dynamic wavelength
    solutions from fitAndApplyEtalonWls. If None, calibration database
    is queried for matching etalon frames. Default is None.

fibers : list of int, optional
    Science fiber numbers to process. Valid values are typically 2, 3, 4
    for MAROON-X science fibers. Default is None.

symmetric_linefits : bool, optional
    If True, use symmetric line profile fitting for etalon peak
    detection. Default is False.

n_knots : int, optional
    Number of interior knots for cubic spline interpolation of
    wavelength solutions. Default is 30.

thar : bool, optional
    Whether to apply ThAr wavelength solution to Etalon frames.
    Default is False. Not implemented.

ref_fiber : int, optional
    Fiber number containing simultaneous etalon spectrum for drift
    measurement. Typically fiber 5 for MAROON-X. Default is None.

report : bool, optional
    Write PDF report with diagnostic plots. Default is True.

suffix : str, optional
    Suffix to append to output filenames. Default is ``'_wls'``.

Returns
-------
list of :class:`~astrodata.AstroData`
    Input frames with the following extension added, one per fiber
    ``N`` in the processed set (science fibers plus ref_fiber):

    - ``WLS_SIMULTANEOUS_FIBER_N`` : 2D array of drift-corrected
      wavelength values (nm) indexed by [order, pixel].

    Header keywords added:

    - ``INSTRUME_DRIFT`` : instrumental drift measured from reference
      fiber in m/s.

    - ``RELATIVE_DRIFT`` : relative drift between science and
      calibration etalon frames in m/s.

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
