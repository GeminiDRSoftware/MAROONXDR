
Perform box extraction on a 2D echelle spectrum.

Utilized in the dynamic and static wavelength calibration recipes
as it is quicker than relying on optimal extraction.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    2D echelle frames carrying the ``STRIPES``, ``F_STRIPES``, and
    ``STRIPES_MASKS`` attributes as dicts of sparse arrays.

suffix : str
    Suffix to be added to output files.

Returns
-------
list of :class:`~astrodata.AstroData`
    The input frames with box-extracted orders for each fiber and
    uncertainties calculated during the box extraction. Extensions
    added, one per fiber ``N`` in ``1-5``:

    - ``REDUCED_ORDERS_FIBER_N``
    - ``BOX_REDUCED_FIBER_N``
    - ``BOX_REDUCED_VAR_N``
    - ``BOX_REDUCED_FLAT_N``
    - ``BPM_FIBER_N``
