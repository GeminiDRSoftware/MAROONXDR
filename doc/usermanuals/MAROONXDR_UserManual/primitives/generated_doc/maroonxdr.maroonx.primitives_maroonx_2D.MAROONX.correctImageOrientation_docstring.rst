
Correct image orientation to proper echelle format.

Flip SCI if needed so that left lower corner is bluest wavelength,
upper right corner is reddest wavelength. Resulting echelle orders
go from left to right. MX blue frames start with incorrect orientation
for reduction. This primitive must be called after DQ is established
and before any image arithmetic is performed.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    List of un-checked MX objects.

suffix : str
    Suffix to be added to output files.

Returns
-------
list of :class:`~astrodata.AstroData`
    Same list as inputs, with correct orientation to SCI.
