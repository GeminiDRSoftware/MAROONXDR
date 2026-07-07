
Check arm consistency across all MAROON-X frames.

Verify that all frames have the same camera arm (BLUE or RED) based on
data tags. The first file sets the expected value. Currently assumes
1 astrodata object comes from 1 single-extension FITS. Need to update
if/when original FITS are MEF.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    List of un-checked MX frames.

suffix : str
    Suffix to be added to output files.

Returns
-------
list of :class:`~astrodata.AstroData`
    Set of frames that pass the test, always at least the first frame.
