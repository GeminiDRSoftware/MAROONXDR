applyBarycentricCorrection
==========================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
| **Astrodata Tags**: {'SCI', 'MAROONX'}

Apply barycentric velocity correction to already-reduced MAROON-X spectra.

Use this recipe to apply target-specific barycentric correction parameters
after the main extraction workflow.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def applyBarycentricCorrection(p):
        p.barycentricCorrection()
        p.storeProcessedScience(suffix='_barycor')

