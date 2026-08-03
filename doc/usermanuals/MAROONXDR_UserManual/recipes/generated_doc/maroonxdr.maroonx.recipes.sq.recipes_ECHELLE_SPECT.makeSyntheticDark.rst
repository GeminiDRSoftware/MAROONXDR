makeSyntheticDark
=================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
| **Astrodata Tags**: {'MAROONX', 'SCI'}

Construct synthetic DDDDE darks for science exposures.

The per-pixel log-linear fit stored in a processed dark coefficients
calibration, retrieved from the calibration database, is evaluated at
the exposure time and ND filter setting of each science frame. This
interpolates the empirical master darks to exposure times that were not
directly observed. The synthetic dark is stored with a "_synth_dark"
suffix.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def makeSyntheticDark(p):
        p.prepare()
        p.checkArm()
        p.addVAR(read_noise=True, poisson_noise=True)
        p.createSyntheticDark()
        p.storeProcessedDark(suffix='_synth_dark')

