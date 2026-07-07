
Separate input AstroData objects into blue and red arm streams.

This primitive takes a list of AstroData objects and separates them
into two lists based on their tags: one for the blue arm and one for
the red arm.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    Input list of AstroData objects to be separated.

Returns
-------
list of :class:`~astrodata.AstroData`
    List of AstroData objects tagged as 'BLUE' arm. The RED arm
    objects are stored in self.streams['RED'].
