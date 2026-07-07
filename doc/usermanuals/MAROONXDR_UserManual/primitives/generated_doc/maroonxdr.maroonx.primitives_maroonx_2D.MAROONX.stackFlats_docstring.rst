
Stack MAROON-X flat field frames.

A simplified DRAGONS primitive that reproduces the legacy make_master_flat.py.

Parameters
----------
adinputs : list of AstroData
    Input frames to be combined
operation : str
    Combine method (default: 'mean')
reject_method : str
    Rejection method (default: 'sigclip')
scale : bool
    Whether to scale images (default: True)

Returns
-------
list of AstroData
    Combined output frame
