makeSyntheticDark
=================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
| **Astrodata Tags**: {'MAROONX', 'SCI'}

Construct DDDDE master darks from coefficient file for science exposures.

Uses a linear interpolation of log(exposure time) vs. flux in empirical
master darks to construct interpolated darks for other exposure times.

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

::

    def makeSyntheticDark(p):
        p.prepare()
        p.checkArm()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.createSyntheticDark()
        p.storeProcessedDark(suffix='_synth_dark')

