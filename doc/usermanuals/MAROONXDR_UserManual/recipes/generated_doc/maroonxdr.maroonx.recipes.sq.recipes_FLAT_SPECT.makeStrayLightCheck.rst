makeStrayLightCheck
===================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
| **Astrodata Tags**: {'CAL', 'FLAT', 'MAROONX'}

Check the stray light subtraction in normal flat frame processing.

Mirrors the makeProcessedFlatDFFFF flow up to removeStrayLight, run with
snapshot=True: the SCI data is left at the level just before straylight
removal and the removed straylight is saved as the STRAYLIGHT_DIFFERENCE
extension. A unit test independently performs the straylight removal and
compares its result against SCI + STRAYLIGHT_DIFFERENCE.

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

    Returns
    -------
    Creates test frames with straylight difference and flux at levels just
    before straylight removal.

::

    def makeStrayLightCheck(p):
        """
        Check the stray light subtraction in normal flat frame processing.

        Mirrors the makeProcessedFlatDFFFF flow up to removeStrayLight, run with
        snapshot=True: the SCI data is left at the level just before straylight
        removal and the removed straylight is saved as the STRAYLIGHT_DIFFERENCE
        extension. A unit test independently performs the straylight removal and
        compares its result against SCI + STRAYLIGHT_DIFFERENCE.

        Parameters
        ----------
        p : PrimitivesCORE object
            A primitive set matching the recipe_tags.

        Returns
        -------
        Creates test frames with straylight difference and flux at levels just
        before straylight removal.
        """
        p.prepare()
        p.checkArm()
        p.addDQ()
        p.subtractOverscan()
        p.addVAR(read_noise=True, poisson_noise=True)
        # Creates 'DFFFD_flats' stream and leaves DDDDF flats in main stream
        p.separateFlatStreams()

        p.stackFlats(stream='main', scale_mode='mean_frame', suffix='DDDDF')
        p.stackFlats(stream='DFFFD_flats', scale_mode='mean_frame', suffix='DFFFD')

        # Subtract overscan is run again in backgroundfit.py on the stacked flats
        p.subtractOverscan(stream='main')
        p.subtractOverscan(stream='DFFFD_flats')

        p.trimOverscan(stream='main')
        p.trimOverscan(stream='DFFFD_flats')

        p.correctImageOrientation(stream='main')
        p.correctImageOrientation(stream='DFFFD_flats')

        # Define stripe info to ultimately remove stray light in each stream
        p.findStripes(stream='main')
        p.findStripes(stream='DFFFD_flats')

        # Identify stripes based on MX architecture files
        p.identifyStripes(stream='main', selected_fibers=[5])
        p.identifyStripes(stream='DFFFD_flats', selected_fibers=[2, 3, 4])

        # Defines pixel inclusion for each flat region based on stripe ids
        p.defineFlatStripes(stream='main')
        p.defineFlatStripes(stream='DFFFD_flats')

        # Remove straylight (requires 2 partial illumination flat sets)
        p.removeStrayLight(stream='main', snapshot=True, filter_size=19, box_size=20)
        p.removeStrayLight(stream='DFFFD_flats', snapshot=True, filter_size=19, box_size=20)

        p.writeOutputs(stream='DFFFD_flats', suffix='_straylight_flat', strip=True)
        p.writeOutputs(suffix='_straylight_flat', strip=True)

