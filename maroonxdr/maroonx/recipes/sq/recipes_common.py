"""
Recipes shared by all MAROON-X recipe libraries.

This module defines no recipe_tags: the recipes here are imported by the
tag-specific recipe libraries, following the GMOS recipes_common example.
"""


def exportReducedBundle(p):
    """
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

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.separateArmStreams()
    p.bundleArmStreams()
    p.writeOutputs()
