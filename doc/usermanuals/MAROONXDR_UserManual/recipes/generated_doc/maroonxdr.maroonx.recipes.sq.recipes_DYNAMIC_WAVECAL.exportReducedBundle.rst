exportReducedBundle
===================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_DYNAMIC_WAVECAL
| **Recipe Imported From**: maroonxdr.maroonx.recipes.sq.recipes_common
| **Astrodata Tags**: {'ARC', 'MAROONX'}

Bundle processed Blue and Red arm frames into a single output file.

Reverses the arm split performed by processBundle: the processed Blue
and Red arm files of the same observation are combined into one
multi-extension bundle, which writeOutputs stores under the name of
the archive file of the observation (e.g. N20250724M4771_reduced.fits)
with a suffix matching the per-arm input products. An explicit suffix
parameter on bundleArmStreams overrides the derived one.

The recipe applies to every processed frame type: reduced science,
processed darks (including dark coefficients and synthetic darks),
processed flats, and dynamic wavelength solutions (etalon or LFC
arcs). ThAr frames are not supported. The recipe bundles one Blue and
Red pair per archive file, so different product types of the same
observation must be exported in separate runs.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

::

    def exportReducedBundle(p):
        p.separateArmStreams()
        p.bundleArmStreams()
        p.writeOutputs()

