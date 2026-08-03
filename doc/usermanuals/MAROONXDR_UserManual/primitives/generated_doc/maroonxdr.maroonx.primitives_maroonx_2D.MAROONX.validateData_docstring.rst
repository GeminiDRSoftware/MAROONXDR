
Validate MAROON-X data while ignoring invalid WCS exceptions.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    List of unchecked AstroData objects.

suffix : str
    Suffix to be added to output files.

require_wcs : bool
    Do all extensions have to have a defined WCS?

Returns
-------
list of :class:`~astrodata.AstroData`
    List of checked AstroData objects.
