
Check ND filter consistency across all frames.

Verify that the ND filter on the sim cal fiber is consistent through
all MX-input files i.e. illumination is similar intensity as needed for
good removal. The first file sets the expected value.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    List of MX-objects.

suffix : str
    Suffix to be added to output files.

Returns
-------
list of :class:`~astrodata.AstroData`
    Inputs that pass the test.
