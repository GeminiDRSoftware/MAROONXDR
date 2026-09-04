
Bundle Blue and Red arm AstroData objects.

This primitive takes the Blue and Red arm streams and combines them
into multi-extension bundle AstroData objects, reversing the operation
performed by splitBundle(). Each bundle contains both Blue and Red arms
as separate extensions.

Parameters
----------
adinputs : list of :class:`~astrodata.AstroData`
    List of Blue arm AstroData objects to be combined with
    previously stored Red arm stream in self.streams['RED'].

suffix : str, optional
    Suffix to append to output filenames. Default is ``'_reduced'``.

Returns
-------
list of :class:`~astrodata.AstroData`
    List of bundle AstroData objects, each containing Blue and Red
    arm extensions with restored ARCHNAME filenames.

Notes
-----
This primitive requires that the separateArmStreams primitive has
been run beforehand to populate self.streams['RED'] with the Red
arm AstroData objects. Each Blue/Red pair must have matching
ARCHNAME headers to be properly bundled together.
