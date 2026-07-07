
Separate flat data into two streams based on fiber setup.

This primitive divides the input flats into two categories:
- 'DFFFD': stored in p.streams['DFFFD_flats']
- 'FDDDF' or 'DDDDF': stored in p.streams['main']

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    MX flats.

suffix : str
    Suffix to be added to output files.

Returns
-------
adoutputs : 'FDDDF' or 'DDDDF' list of AstroData objects

Notes
-----
Modifies the instance's `streams` dictionary by adding a 'DFFFD_flats' key
containing the list of DFFFD flat field AstroData objects.

