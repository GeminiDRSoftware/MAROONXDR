
Extracts the stripes from the original 2D spectrum to a sparse array,
containing only relevant pixels.

This function marks all relevant pixels for extraction.
Reinterpreting the flat reference it iterates over all stripes in the
image and saves a sparse matrix for each stripe.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Input science frames.

suffix : str
    Suffix to be added to output files.

flat : str or :class:`~astrodata.AstroData`, optional
    Processed flat frame. If None, queries the calibration database.

dark : str or :class:`~astrodata.AstroData`, optional
    Processed dark frame. If None, queries the calibration database.

dark_subtraction_skip_fibers : list of int, optional
    Fiber numbers (1-5) to skip dark frame subtraction.

straylight_removal_fibers : list of int, optional
    Fiber numbers (1-5) for which straylight will be removed.

slit_height : int
    Total slit height in px.

test_extraction : bool
    Used in unit test for this function, saves science extraction, flat
    extraction, and the bpm-extraction in FITS-readable format
    (STRIPES, F_STRIPES, STRIPES_MASKS).

report : bool
    Passed along to the straylight removal primitive.

Returns
-------
list of :class:`~astrodata.AstroData`
    Input frames with sparse matrices added holding the 2D extractions
    for each fiber/order for the science frame, flat frame, and BPM.
    Extensions added:

    - ``STRIPES`` : science frame stripes.

    - ``F_STRIPES`` : flat frame stripes.

    - ``STRIPES_MASKS`` : bad pixel mask stripes.

    If ``test_extraction=True``, the extensions are stored in
    FITS-readable format instead of sparse matrix format.
