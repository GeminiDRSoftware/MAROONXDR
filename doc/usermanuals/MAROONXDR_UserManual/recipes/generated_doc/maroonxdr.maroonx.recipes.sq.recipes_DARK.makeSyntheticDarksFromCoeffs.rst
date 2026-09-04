makeSyntheticDarksFromCoeffs
============================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DARK
| **Astrodata Tags**: {'DARK', 'MAROONX', 'CAL'}

Construct synthetic DDDDE darks at given exposure times from coefficients.

The input is a processed dark coefficients calibration (COEFF_Z0 and
COEFF_Z1 extensions) produced by makeDarkCoefficients. For each value
of the exptime list parameter of createSyntheticDarkFromCoeffs, the dark
is evaluated as z1 plus z0 times log10 of the exposure time, and each
result is stored by storeProcessedDark with a "_synth_dark" suffix.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def makeSyntheticDarksFromCoeffs(p):
        p.checkArm()
        p.createSyntheticDarkFromCoeffs()
        p.storeProcessedDark(suffix='_synth_dark')

