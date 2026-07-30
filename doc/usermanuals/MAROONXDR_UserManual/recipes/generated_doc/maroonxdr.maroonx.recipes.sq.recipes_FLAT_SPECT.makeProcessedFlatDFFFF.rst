makeProcessedFlatDFFFF
======================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
| **Astrodata Tags**: {'CAL', 'FLAT', 'MAROONX'}

Perform standardization and corrections to convert raw flats to processed.

This recipe converts the raw input flat images into a single stacked flat
image. This output processed flat is stored on disk using storeProcessedFlat
and has a name equal to the name of the first input bias image with
"_FFFFF_flat.fits" appended. A processed flatfield is required to perform
optimal flux extraction and to determine the blaze function. Due to the
stability of the spectrograph, one processed flatfield is typically 'valid'
for at least a two-week period. Cross-comparison between different processed
flatfields taken months apart have not been conducted so far.

::

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.

::

    def makeProcessedFlatDFFFF(p):
        p.prepare()
        p.checkArm()
        # p.checkND()
        p.addDQ()
        p.subtractOverscan()
        # p.trimOverscan()  # noqa: ERA001
        # p.correctImageOrientation()  # noqa: ERA001
        p.addVAR(read_noise=True, poisson_noise=True)
        # Creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream
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
        # ================================================

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
        p.removeStrayLight(stream='main', filter_size=19, box_size=20)
        p.removeStrayLight(stream='DFFFD_flats', filter_size=19, box_size=20)

        # Legacy patch for removeStrayLight
        # p.removeStrayLight_legacyPatch(stream='main', filter_size=19, box_size=20)
        # p.removeStrayLight_legacyPatch(stream='DFFFD_flats', filter_size=19, box_size=20)

        # Combine straylight-removed images
        p.combineFlatStreams(stream='main', stream_2='DFFFD_flats')
        # Remove second stream
        p.clearStream(stream='DFFFD_flats')

        # Re-run find/identify/define routine on combined frame
        p.findStripes()
        p.identifyStripes(selected_fibers=[2, 3, 4, 5])
        p.defineFlatStripes(extract=True)

        p.storeProcessedFlat(suffix='_DFFFF_flat')

