
Combine two flat field streams into a single unified flat.

This primitive takes two streams of pre-processed flat fields
(typically 'FDDDF' flats in the main stream and 'DFFFD' flats in a
secondary stream) and creates a combined flat field by taking the
maximum value at each pixel position. This produces a flat field with
all fibers illuminated (FFFFF flat configuration).

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Primary stream of flat objects (typically FDDDF flats).

suffix : str
    Suffix to be added to output files.

stream_2 : str
    Name of the secondary stream to combine with the main stream
    (typically 'DFFFD_flats').

Returns
-------
list of :class:`~astrodata.AstroData`
    List containing a single AstroData object with the combined flat
    field data.
