makeStripeExtractionCheck
=========================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
| **Astrodata Tags**: {'MAROONX', 'SCI'}

Check the stripe extraction in normal processing of a science frame.

Test-support recipe. Runs the extraction with the test option enabled so
the stripe extractions, normally kept as sparse arrays, are saved in
FITS-readable format (STRIPES, F_STRIPES and STRIPES_MASKS extensions).
A unit test can then independently perform the stripe extraction and
compare results. The output is written with a "_test_stripes" suffix.

::

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.

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

