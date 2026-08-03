
Stack MAROON-X dark frames with etalon intensity scaling.

MX-specific version of stack darks allowing scaling for etalon
intensity drift that is in MX 'darks'.

Parameters
----------
adinputs : list of AstroData
    Input frames to be combined
suffix : str
    Suffix to be added to output files
scale_mode : str
    Scaling mode for the input frames.
    Options are 'mean_frame' or 'first_frame'.
lsigma : float
    Lower sigma clipping threshold for the rejection method
hsigma : float
    Upper sigma clipping threshold for the rejection method
max_iters : int
    Maximum number of iterations for the rejection method
reject_method : str
    Rejection method to be used. Currently only 'sigclip' is supported.

Returns
-------
list of AstroData
    Combined output frame
