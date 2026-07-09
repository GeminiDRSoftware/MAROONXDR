makeStrayLightCheck
===================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
| **Astrodata Tags**: {'FLAT', 'MAROONX', 'CAL'}

Check the stray light subtraction in normal flat frame processing.

Run the straylight_test_prep.py file to generate these.

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

    Returns
    -------
    Creates test frames with straylight difference and flux at levels just
    before straylight removal. Unit test will independently perform straylight
    removal and compare results.

::

    def makeStrayLightCheck(p):
        p.prepare()
        p.checkArm()
        # p.checkND()
        p.addDQ()
        p.subtractOverscan()
        p.trimOverscan()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        # Creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
        p.separateFlatStreams()
        p.stackFlats(suffix='FDDDF_flats')
        p.stackFlats(stream='DFFFD_flats', suffix='DFFFD')
        # Define stripe info to ultimately remove stray light in each stream
        p.findStripes()
        p.findStripes(stream='DFFFD_flats')
        # Identify stripes based on MX architecture files
        p.identifyStripes(selected_fibers=[1, 5])
        p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])
        # Defines pixel inclusion for each flat region based on stripe ids
        p.defineFlatStripes()
        p.defineFlatStripes(stream='DFFFD_flats')
        # Remove straylight (requires 2 partial illumination flat sets)
        p.removeStrayLight(snapshot=True, filter_size=19, box_size=20)
        p.removeStrayLight(stream='DFFFD_flats', snapshot=True, filter_size=19, box_size=20)
        p.writeOutputs(stream='DFFFD_flats', suffix='_straylight_flat')
        p.writeOutputs(suffix='_straylight_flat')

