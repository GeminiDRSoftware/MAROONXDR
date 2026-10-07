
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
    Suffix appended to output filenames. If None (the default),
    the suffix is derived per bundle from the bundle tags so the
    output name mirrors the per-arm input products, checking the
    most specific tag first: ``DARK_COEFF`` gives
    ``'_darkCoefficients'``, ``DARK_SYNTH`` gives ``'_synth_dark'``,
    ``DARK`` gives ``'_dark'``, ``FLAT`` gives ``'_flat'``, and
    ``ARC`` gives ``'_arc'``; ``'_reduced'`` is used for science
    frames and as the fallback when no mapped tag is present.
    An explicitly passed suffix, including an empty string,
    always overrides the derivation.

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
ARCHNAME headers to be properly bundled together, and exactly one
bundle is produced per ARCHNAME. If several files of the same arm
share an ARCHNAME (e.g. a processed dark and its dark
coefficients), only the first one in the input list is kept and a
warning is logged for each skipped file. To bundle mixed product types, feed the
primitive one product type at a time.
