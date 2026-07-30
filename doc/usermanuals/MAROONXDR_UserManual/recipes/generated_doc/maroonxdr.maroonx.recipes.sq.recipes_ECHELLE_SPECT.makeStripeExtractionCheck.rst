makeStripeExtractionCheck
=========================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
| **Astrodata Tags**: {'MAROONX', 'SCI'}

Check the stripe extraction in normal processing of a science frame.

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

    Returns
    -------
    Creates test frames with FITS-formatted stripe extractions meta-info
    (normally not saved). Unit test will independently perform stripe extraction
    and compare results.

::

    def makeStripeExtractionCheck(p):
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        # Gets relevant flat and dark to cut out frame's spectra
        p.extractStripes(
            dark_subtraction_skip_fibers=[5],
            straylight_removal_fibers=[5],
            test_extraction=True,
        )
        p.writeOutputs(suffix='_test_stripes')

