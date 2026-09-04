
Create synthetic dark frames matching the exposure time of the inputs.

The per-pixel log-linear model of a processed dark coefficients
calibration, dark = z1 + z0 * log10(exptime), is evaluated at the
exposure time of each input. By default inputs are grouped by
exposure time and arm, and one dark is produced per group.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Frames to create synthetic darks for. Exposure time, ND filter
    position and arm tag (BLUE or RED) are read from each frame.

dark_coeff : str or :class:`~astrodata.AstroData`
    Processed dark coefficients calibration. If None, it is retrieved
    from the calibration database for each input.

individual : bool
    If False, one synthetic dark per unique exposure time and arm,
    named after the first frame of each group. If True, one synthetic
    dark per input frame.

Returns
-------
list of :class:`~astrodata.AstroData`
    Synthetic dark frames only; the inputs are not returned. Each is a
    deep copy of the first frame of its group with the data replaced
    by the evaluated dark and the PHU fiber keywords set to
    ``FIBER1`` to ``FIBER4`` = ``'Dark'``, ``FIBER5`` = ``'Etalon'``.
    Inputs already covered by a previous group, or with no dark
    coefficients available, produce no output.
