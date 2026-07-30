exportReducedBundle
===================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
| **Astrodata Tags**: {'MAROONX', 'SCI'}

Export reduced MAROON-X spectra from Red and Blue channels to bundle.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def exportReducedBundle(p):
        p.separateArmStreams()
        p.bundleArmStreams()
        p.storeProcessedScience(suffix='_reduced')

