makeProcessedFlatDFFFF
======================

| **Recipe Library**: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
| **Astrodata Tags**: {'CAL', 'MAROONX', 'FLAT'}

Convert raw MAROON-X flat frames into a processed DFFFF flat.

Variant of makeProcessedFlat for datasets without FDDDF flats: the
DDDDF flats (only fiber 5 illuminated) are kept in the main stream and
the DFFFD flats (fibers 2, 3 and 4 illuminated) in a second stream.
Following the legacy processing order, each stream is stacked first;
the stacked frames are then overscan subtracted again, trimmed and
orientation corrected. The fiber stripes of each stream are traced
and identified, its stray light is removed, and the two streams are
combined into a DFFFF flat (fiber 1 dark) by taking the pixel-by-pixel
maximum. The stripe tracing is re-run on the combined frame and the
result is stored on disk by storeProcessedFlat with a "_DFFFF_flat"
suffix.

::

    Parameters
    ----------
    p : Primitives object
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

