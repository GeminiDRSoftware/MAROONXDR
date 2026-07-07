
Extract etalon peak positions and fit line profile parameters.

TODO: Rewrite docstring to make reference to maroonx_fit module

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Input AstroData objects with 1D box extracted etalon spectra.
    Must have ETALON tag.

guess_file : str, optional
    Path to file containing initial guess spectrum for peak positions.
    Default is None.

fibers : list of int, optional
    Fiber numbers to process. Valid values are 1-5. If None,
    automatically detects fibers with ETALON or LFC fiber type.
    Default is None.

orders : list of int, optional
    Order numbers to process. If None, all orders are processed.
    Default is None.

degree_sigma : int, optional
    Polynomial degree for Gaussian sigma variation along order.
    Default is 4.

degree_width : int, optional
    Polynomial degree for box width variation along order.
    Default is 2.

use_sigma_lr : bool, optional
    If True, use different polynomial coefficients for left and right
    Gaussian wings. Default is True.

show_plots : bool, optional
    If True, generate diagnostic plots of etalon line fits.
    Automatically disabled if multithreading is True. Default is False.

plot_path : str, optional
    Directory path for saving diagnostic plots when show_plots is True.
    Default is empty string (current directory).

multithreading : bool, optional
    If True, use multiprocessing to parallelize fiber/order fitting.
    Disables show_plots option. Default is False.

iterations : int, optional
    Maximum number of iterative fitting cycles. Default is 8.

suffix : str, optional
    Suffix to append to output filenames. Default is empty string.

Returns
-------
list of :class:`~astrodata.AstroData`
    Input frames with the following extensions added:

    - ``PEAKS`` : Peak parameters including centroid positions,
      intensities, and widths for each detected etalon line.

    - ``POLY`` : Polynomial coefficients describing how box width
      and Gaussian sigma vary across each order.

Notes
-----
The iterative fitting algorithm may fail for orders with excessive
cosmic rays or artifacts. Failed fits are logged as warnings.

Known bad pixels are masked during fitting (e.g., pixel 1943 in order
122, first 400 pixels in truncated red arm order 94 fiber 5).
