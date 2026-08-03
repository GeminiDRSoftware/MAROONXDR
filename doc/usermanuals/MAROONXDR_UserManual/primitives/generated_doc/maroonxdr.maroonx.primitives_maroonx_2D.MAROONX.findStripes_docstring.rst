
Locate and fit stripes in flat field spectrum.

Starting in the central column, the algorithm identifies peaks and
traces each stripe to the edge of the detector by following the
brightest pixels along each order. It then fits a polynomial to each
stripe. To improve algorithm stability, the image is first median
filtered and then smoothed with a gaussian. It not only eliminates
noise, but also ensures that the cross-section profile of the flat
becomes peaked in the middle, which helps to identify the center of
each stripe. Choose ``gauss_filter_sigma`` accordingly. To avoid
false positives, only peaks above a certain (relative) intensity
threshold are used.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    MX flat frames. Each frame is either a DFFFD flat, an FDDDF
    flat, or a combined FFFFF flat.

suffix : str
    Suffix to be added to output files.

deg_polynomial : int
    Degree of the polynomial fit to each stripe.

med_filter : int
    Median filter window size applied before peak detection.

gauss_filter_sigma : float
    Sigma of the gaussian filter used to smooth the image before
    peak detection.

min_peak : float
    Minimum peak height, relative to the frame maximum, for a
    stripe to be accepted.

Returns
-------
list of :class:`~astrodata.AstroData`
    The input frame with the following extension added:

    - ``STRIPES_LOC`` : per-fiber polynomial coefficients tracing
      each stripe. This extension temporarily holds the
      fits-unsavable fiber information before it is utilized and
      then removed.
