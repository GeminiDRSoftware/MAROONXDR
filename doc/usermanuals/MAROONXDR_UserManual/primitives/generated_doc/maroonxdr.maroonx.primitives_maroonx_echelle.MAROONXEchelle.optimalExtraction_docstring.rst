
Optimal extraction of the 2d echelle spectrum.

This function performs an optimal extraction of a 2d echelle spectrum.
The internal flat spectrum extension is used to generate normalized
'profiles' that are used as weighting functions for the spectrum that
is going to be extracted.
The algorithm further checks for outliers and rejects them.
This is to prevent contributions from cosmic hits.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Frames with STRIPES, F_STRIPES, and STRIPES_MASKS extensions as
    dicts of sparse arrays.

suffix : str
    Suffix to be added to output files.

dark : str or :class:`~astrodata.AstroData`, optional
    Processed dark frame. If None, queries the calibration database.

optimal_extraction_fibers : list of int, optional
    Fiber numbers (1-5) for optimal extraction. Fibers considered for
    optimal extraction.

back_var : float, optional
    Manual background variance for frame.

full_output : bool
    If True, an additional set of intermediate products will be
    returned / saved.

penalty : float
    Scaling penalty factor for mismatch correction between flat field
    profile and science spectrum during optimal extraction.

s_clip : float
    Sigma-clipping parameter during optimal extraction.

read_noise : float
    Read noise (e-). MX default is 1.14.

gain : float
    Detector gain (e-/ADU). MX default is 2.72.

Returns
-------
list of :class:`~astrodata.AstroData`
    Input frames with optimal and box extracted orders for each fiber
    together with uncertainties and the bad pixel mask result from the
    optimal extraction. Extensions added, one per fiber ``N`` in ``1-5``:

    - ``REDUCED_ORDERS_FIBER_N``
    - ``OPTIMAL_REDUCED_FIBER_N``
    - ``OPTIMAL_REDUCED_VAR_N``
    - ``BOX_REDUCED_FIBER_N``
    - ``BOX_REDUCED_VAR_N``
    - ``BOX_REDUCED_FLAT_N``
    - ``BPM_FIBER_N``
