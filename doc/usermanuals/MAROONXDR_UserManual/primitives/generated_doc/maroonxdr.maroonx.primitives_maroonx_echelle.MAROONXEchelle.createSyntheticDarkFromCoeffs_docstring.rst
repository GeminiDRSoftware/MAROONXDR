
Create synthetic dark frames at given exposure times from coefficients.

The per-pixel log-linear model of the input processed dark
coefficients calibration, dark = z1 + z0 * log10(exptime), is
evaluated at each requested exposure time. Each product is a copy of
the input with the data replaced, the COEFF_Z0, COEFF_Z1 and
LOGEXPTIME extensions removed, EXPTIME set in the primary and
extension headers, and the exposure time field of the filename
rewritten (also in ORIGNAME).

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Processed dark coefficients calibrations (tag DARK_COEFF).

exptime : float or list of float
    Exposure times in seconds, one synthetic dark per value.

Returns
-------
list of :class:`~astrodata.AstroData`
    Synthetic dark frames, one per input and exposure time.

Raises
------
ValueError
    If ``exptime`` is empty or an input is not tagged DARK_COEFF.
