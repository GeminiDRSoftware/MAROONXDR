applyBarycentricCorrection
==========================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
| **Astrodata Tags**: {'SCI', 'MAROONX'}

Apply barycentric velocity correction to already reduced MAROON-X spectra.

Use this recipe to recompute the barycentric correction with target
specific parameters (SIMBAD name, telescope coordinates, exposure meter
zeropoints) after the main extraction workflow. The computed BERV values
and timing information are stored as header keywords and the result is
written with a "_barycor" suffix.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def applyBarycentricCorrection(p):
        p.barycentricCorrection()
        p.storeProcessedScience(suffix='_barycor')

