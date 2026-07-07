
Load static wavelength solution from lookup file.

This primitive loads pre-computed static wavelength calibration solutions
from lookup files and attaches them as an extensions.
The static solutions provide initial wavelength mappings before dynamic
refinement using etalon measurements.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Input AstroData objects with 1D box extracted spectra. Must have
    REDUCED_ORDERS_FIBER_* extensions for each fiber.

fibers : list of int, optional
    Fiber numbers to load wavelength solutions for. If None, all
    fibers (1-5) are processed. Default is None.

suffix : str, optional
    Suffix to append to output filenames. Default is empty string.

Returns
-------
list of :class:`~astrodata.AstroData`
    Input frames with the following extension added, one per fiber
    ``N`` in ``1-5``:

    - ``WLS_STATIC_FIBER_N`` : 2D array of wavelength values (nm)
      indexed by [order, pixel].
