
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
adinputs : list of :class:`~astrodata.AstroData`
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

report : bool, optional
    Write PDF report with diagnostic plots. Default is True.

suffix : str, optional
    Suffix to append to output filenames. Default is ``'_etalonwls'``.

Returns
-------
list of :class:`~astrodata.AstroData`
    Input frames with the following extensions added, one per fiber
    ``N`` in ``1-5``:

    - ``WLS_DYNAMIC_FIBER_N`` : 2D array of wavelength values (nm)
      indexed by [order, pixel].

    - ``PEAK_DATA`` : updated table with wavelength assignments and
      peak order numbers.

    Header keywords added, one per fiber ``N`` in ``1-5``:

    - ``DRIFT_FIBER_N`` : measured instrumental drift in m/s.

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
