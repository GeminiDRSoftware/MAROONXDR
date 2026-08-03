
Fit a blaze function to each fiber of a processed masterflat and store
the results as named AstroData extensions.

Instantiates a :class:`FlatSpectrum` per fiber from the
``BOX_REDUCED_FLAT_{f}`` and ``REDUCED_ORDERS_FIBER_{f}`` extensions
produced by :meth:`optimalExtraction`, calls
:meth:`~FlatSpectrum.fit_blaze`, then stores the resulting normalised
blaze as ``BLAZE_FIBER_{f}`` — a 2-D float32 array of shape
``[n_orders, n_pixels_per_order]`` with ``max == 1`` per order row.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Processed flat AstroData objects that contain
    ``BOX_REDUCED_FLAT_{f}`` extensions.

suffix : str
    Filename suffix appended to output files.

n_knots : int
    Number of spline knots used in the blaze fit (default 50).

outlier_threshold : float
    Relative outlier rejection threshold (fraction of blaze fit value).

fibers : list of int
    Fiber numbers to process.

Returns
-------
list of :class:`~astrodata.AstroData`
    Input frames with the following extension added, one per fiber
    ``N`` in ``1-5``:

    - ``BLAZE_FIBER_N``
