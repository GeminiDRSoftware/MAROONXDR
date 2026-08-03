
Identify stripes by assigning order and fiber numbers.

Assign proper order and fiber number to each stripe, including
correction for the possibility that the spectra have shifted up/down
in the cross-dispersion direction since the reference was made.
Requires the findStripes primitive to be run prior during recipe so
the stripes are located in the input, i.e. STRIPES_LOC extension exists.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    MX flat frames with STRIPES_LOC extension. Each frame is
    either a DFFFD flat, an FDDDF flat, or a combined FFFFF flat.

suffix : str
    Suffix to be added to output files.

positions_dir : str
    Lookup fits location of nominal y positions and fiber/order
    labels. Shape is Nx3, columns are [fibers, orders, y_positions],
    nominally found in lookups/SID.

selected_fibers : list of int
    List of fiber numbers illuminated in the flat. If None, assumes
    all. Can work if not given on partially illuminated frame, but
    best practice is to explicitly identify on function call.

Returns
-------
list of :class:`~astrodata.AstroData`
    The input frames with the following extensions added:

    - ``STRIPES_ID`` : per-fiber, per-order stripe identification.
      This extension temporarily holds the fits-unsavable fiber
      information before it is utilized and then removed.

    - ``REMOVED_STRIPES`` : polynomial info for every stripe that
      is not identified, from the original set inherited from
      findStripes.
