
Define flat field stripe locations for extraction and straylight removal.

Save fiber location info based on flat field info for stray light
removal (extract=False) and for future science extraction (extract=True).
Requires the findStripes and identifyStripes primitives to be run prior
during recipe so necessary information exists in input extensions.
Will remove previous (improperly formatted, but fast)
STRIPES_ID and STRIPES_LOC extensions and replace with INDEX_FIBER
and INDEX_ORDER pixel map extensions, as needed in straylight removal,
and (if extract=True) a FITS savable STRIPES_ID and STRIPES_FIBERS.

For a given slit_height, this function extracts the flat field stripes,
calculates a box extracted spectrum and normalizes the flat field to
generate a 2D pixel map that is used in the straylight removal.

STRIPES_ID and STRIPES_FIBERS contain the by-spectral-order
polynomial plate solution for each illuminated fiber that is utilized
to define 2D extraction regions in science extractions.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    MX flat frames. Each frame is either a DFFFD flat, an FDDDF flat,
    or a combined FFFFF flat.

suffix : str
    Suffix to be added to output files.

slit_height : int
    Half pixel height of box in each dimension to perform box
    extraction with.

extract : bool
    If True, will write STRIPES_ID in fits-acceptable format. Utilized
    in combined, all fiber illuminated FFFFF_flat.

Returns
-------
list of :class:`~astrodata.AstroData`
    The input frames with the following extensions added:

    - ``INDEX_FIBER`` : 2D pixel map of fiber assignments.

    - ``INDEX_ORDER`` : 2D pixel map of echelle order assignments.

    - ``STRIPES_ID`` : (if ``extract=True``) FITS-savable
      by-spectral-order polynomial plate solution for each
      illuminated fiber.

    - ``STRIPES_FIBERS`` : (if ``extract=True``) companion to
      ``STRIPES_ID`` used to define 2D extraction regions in
      science extractions.
