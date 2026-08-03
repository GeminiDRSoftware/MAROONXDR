
Check that MX frames are processed master frames.

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
