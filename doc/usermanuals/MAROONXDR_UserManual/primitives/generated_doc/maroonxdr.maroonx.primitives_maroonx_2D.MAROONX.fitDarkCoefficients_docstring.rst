
Construct coefficients for log-linear fit of flux vs. exposure time.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Input frames to be combined.

suffix : str
    Suffix to be added to output files.

Returns
-------
list of :class:`~astrodata.AstroData`
    Single-element list containing the combined output frame with the
    following extensions added:

    - ``COEFF_Z0`` : 2D array of the log-linear slope coefficients.

    - ``COEFF_Z1`` : 2D array of the log-linear intercept coefficients.

    - ``LOGEXPTIME`` : Table with columns ``logexptime``, ``exptime``,
      ``filename``.
