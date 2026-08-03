
Add variance extension to MAROON-X frames.

Calculate the variance based on the read noise for the chip and the
poisson noise (the variance in this case is just the number of photons
for each pixel). The variance is then stored as a FITS extension for
each file.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    List of MX objects without variance extensions.

suffix : str
    Suffix to be added to output files.

read_noise : bool
    Whether to include read noise in variance calculations.

poisson_noise : bool
    Whether to include poisson noise in variance calculations.

Returns
-------
list of :class:`~astrodata.AstroData`
    List of MX objects with variance extensions.
