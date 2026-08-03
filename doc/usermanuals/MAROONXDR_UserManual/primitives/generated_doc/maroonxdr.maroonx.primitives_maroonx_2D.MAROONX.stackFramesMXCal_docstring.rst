
Stack MAROON-X calibration frames with etalon flux scaling.

MX-specific version of stackFrames for calibration frames - changes
scaling to average full frame mean to purposely scale by etalon flux
and its drift between calibration exposures. This function should
not be used to combine MX science frames.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Any set of 2D.

suffix : str
    Suffix to be added to output files.

apply_dq : bool
    Apply DQ mask to data before combining?

nlow, nhigh : int
    Number of low and high pixels to reject, for the 'minmax' method.
    The way it works is inherited from IRAF: the fraction is specified
    as the number of  high  and low  pixels,  the  nhigh and nlow
    parameters, when data from all the input images are used.  If
    pixels  have  been  rejected  by offseting,  masking, or
    thresholding then a matching fraction of the remaining pixels,
    truncated to an integer, are used.  Thus::

        nl = n * nlow/nimages + 0.001
        nh = n * nhigh/nimages + 0.001

    where n is the number of pixels  surviving  offseting,  masking,
    and  thresholding,  nimages  is the number of input images, nlow
    and nhigh are task parameters  and  nl  and  nh  are  the  final
    number  of  low  and high pixels rejected by the algorithm.  The
    factor of 0.001 is to adjust for rounding of the ratio.

operation : str
    Combine method.

reject_method : str
    Pixel rejection method (none, minmax, sigclip, varclip).

zero : bool
    Apply zero-level offset to match background levels?

scale : bool
    Scale images to the same intensity?

memory : float or None
    Available memory (in GB) for stacking calculations.

statsec : str
    Section for statistics.

separate_ext : bool
    Handle extensions separately?

Returns
-------
list of :class:`~astrodata.AstroData`
    Sky stacked image. This list contains only one element. The list
    format is maintained so this primitive is consistent with all the
    others.

Raises
------
IOError
    If the number of extensions in any of the `AstroData` objects is
    different.

IOError
    If the shape of any extension in any `AstroData` object is different.

AssertError
    If any of the `.gain()` descriptors is None.

AssertError
    If any of the `.read_noise()` descriptors is None.
