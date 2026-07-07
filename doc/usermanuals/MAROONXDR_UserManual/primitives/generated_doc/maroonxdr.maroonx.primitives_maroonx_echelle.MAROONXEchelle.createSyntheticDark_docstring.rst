
Creates synthetic dark frames from science input files.

This primitive generates synthetic dark frames using pre-computed coefficients
and the log-linear relationship: dark = z1 + z0 * log10(exptime * factor)

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Science files to create darks for.

suffix : str
    Suffix to be added to output files.

dark_coeff : str or :class:`~astrodata.AstroData`, optional
    Adinput of dark coefficients file.

individual : bool
    If True, creates unique dark for each frame (no reuse). If False,
    reuses darks for frames with same exposure time, ND filter, and arm.

Returns
-------
list of :class:`~astrodata.AstroData`
    Synthetic dark frames.
