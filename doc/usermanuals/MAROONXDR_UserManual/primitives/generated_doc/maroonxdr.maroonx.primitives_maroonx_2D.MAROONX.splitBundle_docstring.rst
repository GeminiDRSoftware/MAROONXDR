
Split bundle into separate AstroData objects per extension.

Split a bundle (Red and Blue multi-extension AstroData object) into
separate AstroData objects, each containing one of the original extensions.

This primitive creates a separate AstroData object for each extension
in the input bundle. Each new object is a complete deep copy of the
original, with all other extensions removed. This ensures all tables
and associated data are preserved for the remaining extension. The
original file name is retrieved from the 'ORIGNAME' header and assigned
to the new object, and an 'ARCHNAME' header is added to reference the
Gemini Archive file name. Additionally, certain header keywords in the
EXPOSUREMETER table are renamed to avoid being lost by writing methods.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    List of multi-extension AstroData objects to be split.

suffix : str
    Suffix to be added to output files.

keep_suffix : bool
    If True, retains the suffix from the processed file name.

Returns
-------
list of :class:`~astrodata.AstroData`
    List containing a separate AstroData object for each extension in
    each input.
