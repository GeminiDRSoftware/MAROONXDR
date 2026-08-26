exportReducedBundle
===================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
| **Astrodata Tags**: {'MAROONX', 'SCI'}

Bundle reduced Red and Blue arm spectra into a single output file.

Reverses the arm split performed by processBundle: the reduced Blue and
Red arm files of the same observation are combined into one
multi-extension bundle, which is stored with a "_reduced" suffix.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def exportReducedBundle(p):
        p.separateArmStreams()
        p.bundleArmStreams()
        p.storeProcessedScience()

