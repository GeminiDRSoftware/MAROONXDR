# MAROONXDR

This repo contains the MAROON-X DRAGONS implementation of the data reduction pipeline. We recommend reading this document in its entirety. **PLEASE ENSURE YOU AT LEAST READ THE QUICK INSTALL SECTION BEFORE REPORTING ANY ERRORS!**

## TABLE OF CONTENTS

- [QUICK INSTALLATION](#quick-installation)
- [INTRODUCTION](#introduction)
- [EXPECTATIONS FOR THE RAW FRAMES PROVIDED TO THE PIPELINE](#expectations-for-the-raw-frames-provided-to-the-pipeline)
- [RUN OPTIONS](#run-options)
- [RECIPE DETAILS](#recipe-details)
  - [MASTER DARK CREATION](#master-dark-creation)
  - [MASTER FLAT CREATION](#master-flat-creation)
  - [SCIENCE FLUX EXTRACTION](#science-flux-extraction)
- [DETAILS ON PRIMITIVES IN PRIMITIVES_MAROONX.PY](#details-on-primitives-in-primitives_maroonxpy)
- [DETAILS ON PRIMITIVES IN PRIMITIVES_MAROONX_ECHELLE.PY](#details-on-primitives-in-primitives_maroonx_echellepy)
- [IMAGE TESTING](#image-testing)

## QUICK INSTALLATION

There are certain non-obvious steps that are currently necessary to ensure Dragons and Maroon-X data reduction work properly.
Not following these steps will probably lead to run-time errors that are very difficult to diagnose.

1. If not already installed, install Anaconda.  If on an M1 Mac, ensure that you are using the x86 version of the Anaconda binaries.  The M1 version of Anaconda does not work with DRAGONS.  The x86 version can be found here: <https://www.anaconda.com/products/individual>

2. Set up a new conda environment with python version 3.10, Dragons, DS9, astropy version 5.2 and bokeh version 2.4.3.   The command will be formatted like this: conda create -n NAME python=3.10 dragons ds9 astropy=5.2 bokeh=2.4.3 <br>
**NOTE:  IT IS CURRENTLY ESSENTIAL TO ENSURE YOU INSTALL OLDER VERSIONS OF ASTROPY AND BOKEH.  DRAGONS DOES NOT WORK WITH THE LATEST VERSION OF ASTROPY.**

3. Edit dark_test.py to provide the correct path for your dark files (DDDDE files), and run this file.  It should execute correctly.  Following this, do the same for flat_test.py and science_test.py.  If all three files execute correctly, you have installed DRAGONS correctly.

The simplest way to ensure correct installation is to create a science_dir file in the root folder of the installation and put all fits files in this folder.  Doing this will ensure that you do not have to change the path in simple_test.py.

## INTRODUCTION

MAROON-X is a fiber-fed echelle-spectrograph installed on the Gemini-North Observatory meant to detect Earth-sized planets in the habitable zones of mid- to late-M dwarves by measuring their radial velocities with 1 m/s radial velocity precision.

The red-optical, fiber-fed echelle-spectrograph has no movable parts and only one mode. At each observed wavelength the 5 fibers are arranged along the cross-dispersion direction. The middle three fibers are pupil-sliced fractions of the on-sky target fiber, the outer two fiber traces are an off-target on-sky fiber and a fiber known as the 'sim. cal. fiber' which can be fed calibration light during science frames (as well as during other frames). All observed light is split with a dichroic to a ‘blue arm’ (491-670nm) and a ‘red arm’ (649-920nm). Both arms terminate in CCD detectors with 16-bit 4400x4400 arrays, but the 'blue' detector has 4 reads and the 'red' detector has 2.

The repo is currently a re-implementation of the 'flux extraction' part of the internal U-Chicago pipeline. This first step occurs before a second wavelength calibration step of the internal two-step pipeline.  There are 2 types of wavelength calibrations- dynamic and static.  Dynamic wavelength calibrations are carried out every night of observing, while static wavelength calibrations have only been carried out once in the instrument's lifetime.  For this reason, this repo will probably not include recipes for static wavelength calibration- this remaining only in the internal U-Chicago pipeline.

Broadly the repo currently makes a few products that are recognizable as standard products:

1. Master Darks
2. Master Flats
3. 2D echelle optimally extracted reduced science frames

In progress are the wavelength calibration recipes : DYNAMICWAVECAL, the end of the echelle_spect recipe, and the to be written additional primitive calls in the echelle_spect recipe.

## EXPECTATIONS FOR THE RAW FRAMES PROVIDED TO THE PIPELINE

Repo currently expects on any raw MAROON-X fits read:

1. All raw maroon-x fits files to have keyword INSTRUME = 'MAROON-X' (using the matches_data method repo then internally removes hyphen in instrument definition to agree with repo names)

2. All raw maroon-x fits files with corresponding red and blue arm exposures as completely separate files (currently based on file name in data_tag: _tag_arm)

3. All raw fits frames are full-frame images with no changing bias/data/array/detector section values. Each read section is treated as a separate data section with a unique name. The values are stored in the ./maroonx_instruments/maroonx/lookups.py file.

4. All raw maroon-x fits files to have defined HIERARCH FIBER1 through FIBER5 for all frames (multiple data_tags utilize these entries to to set the associated TagSet for any file type)
-- currently utilized data_tags are _tag_dark,_tag_flat, and _tag_science which currently utilize 'Dark'+'Etalon', 'Flat', and 'OBJECT'+'ETALON' entry values respectively to make 'DARK'+'CAL', 'FLAT'+'CAL', and 'SCI' TagSet associations respectively.

5. Currently with any raw frame, HIERARCH MAROONX IMAGE ORIENTATION HORIZONTAL FLIP and HIERARCH MAROONX IMAGE ORIENTATION VERTICAL FLIP are utilized to re-orient blue frames to have the same dispersion and cross-dispersion that red frames have. This is redundant with _tag_arm, and we will make the change soon to use that instead.

## RUN OPTIONS

Performing a standard 'reduce' command on all frames of a particular arm and file type can be achieved through the standard methods.  Please refer simple_test.py to see an example of how a reduction can be run.  
DRAGONS requires the user to ensure that the correct files are being provided- it will only check the tags of the first image to see what recipe should be used, and any frame with tags that do not match the initial frame will simply be tossed.  Use the select_data() function provided by the DRAGONS API to ensure you do not provide incorrect frames.

## RECIPE DETAILS

### MASTER DARK CREATION

- found at ./maroonxdr/maroonx/sq/recipes_DARK.py

This recipe performs the standardization and corrections needed to convert the raw input dark images into a single processed dark image. This output processed dark is stored on disk using storeProcessedDark and has a name equal to the name of the first input bias image with "_dark.fits" appended. The background in an un-illuminated frame is very low for exposure times of less than 900s and likely doesn't warrant a dark subtraction, however most science frames are taken with the simultaneous calibration fiber illuminated with the FP etalon. The extended wings of the etalon reach into one of the science fibers with a few 10 counts. To remove these wings and any broad diffuse (illumination independent) background, DDDDE frames are taken in daytime to construct a DDDDE processed dark. These darks are specific for different exposure times (i.e. ND filter settings) and should be taken close in time (within a day or two) to the science frame as the etalon source brightness can be time variable.

Notes on recipe in chronological order of primitive calls.  Please refer to the list of primitives for a more detailed explanation about what each primitive does.

- checkArm: Robust to being fed red and blue arm frames.  Will only use the frames in agreement with the first dark's arm and will log a warning if frames with the incorrect arm is given.  If only one frame is given, an error will be raised.

- checkND: Robust to being fed flats at different exposure times with the checkND primitive:
will only use the flats in agreement with the first dark's ND and will log warning, error if only one

- addDQ: Adds BPM in ./maroonxdr/maroonx/lookups/BPM/ depending on data series arm with the addDQ overwrite of standard call in maroonx primitives (only needed because we are not in caldb). Differences are noted with !! commments.  This will change once MAROON-X is fully integrated with DRAGONS.

- overscanCorrect: Overscan corrects each frames depending on arm and instrument lookup values.

- correctImageOrientation: Corrects image orientation, as needed with blue arm, so that left lower corner is bluest wavelength, upper right corner is reddest wavelength. As implemented, this goes (as needed) after overscan correction because it changes blue arm pixel index identification. Refer to the [Details](#details-on-primitives-in-primitives_maroonxpy) section of this file for more information.

- stackDarks: Stacks darks using maroonx primitive overwrite of standard. Refer to the [Details](#details-on-primitives-in-primitives_maroonxpy) section of this file for more information.

- storeProcessedDark: Standard reduction on the recipe spits out the final frame in ./calibrations/processed_dark/

### MASTER FLAT CREATION

- found in default at ./maroonxdr/maroonx/sq/recipes_FLAT_SPECT.py
(for makeStrayLightCheck info, see test_stray_light_removal.py test)

This recipe performs the standardization and corrections needed to convert the raw input flat images into a single stacked flat image. This output processed flat is stored on disk using storeProcessedFlat and has a name equal to the name of the first input bias image with "_FFFFF_flat.fits" appended. A processed flatfield is required to perform optimal flux extraction and to determine the blaze function. Due to the stability of the spectrograph, one processed flatfield is typically 'valid' for at least a two-week period. Cross-comparison between different processed flatfields taken months apart have not been conducted so far. The recipe expects two sets of sub-illuminated flat frames (as needed for best background subtraction processing).

notes on recipe in chronological order of primitive calls:

Performs checkArm, checkND, addDQ, overscanCorrect, and correctImageOrientation identically to their use in the creation of a master dark to ensure all received raw flat frames have identical setup and are prepared correctly for use.

- separateFlatStreams: uses fiber_setup data_descriptor to identify which of the given frames have flat illumination in the non-target fibers (known internally in the code as FDDDF) and stay in the main stream and which have flat illumination in the target fibers (i.e. DFFFD) that get moved to the newly created DFFFD_flats stream. Warns if no frames are found in either of the two patterns.

**Starting here primitive calls are run on both the main stream and the 'DFFFD_flats' stream**

- stackFlats: Stacks based on the stackFramesMXCal primitive that scales based on average full-frame mean.

- findStripes: Locates and fits stripes in the stacked flat field spectrum.  See [Details](#details-on-primitives-in-primitives_maroonxpy) section of this file for more information.

- identifyStripes: Identifies the stripes to their proper order and fiber number, including correction for the possibility that the spectra have shifted up/down in the cross-dispersion direction since the reference was made.  See [Details](#details-on-primitives-in-primitives_maroonxpy) section of this file for more information.

- defineFlatStripes: The first two times this primitive is called (extract==Default==False) is to map all pixels in each stacked flat field frame that have been identified as part of a real fiber trace so that they can be masked in the following background modelling. Requires the findStripes and identifyStripes primitives to be run prior during recipe so necessary information exists in input extensions.  Previous, temporary non-fits-writeable extensions are dropped. I.E. with extract == False creates two 2D composite map frames with the pixel entries being their fibers and orders respectively for all pixels found (by identifyStripes) to be in a fiber trace.  See [Details](#details-on-primitives-in-primitives_maroonxpy) section of this file for more information.

- removeStrayLight: Removes stray light from full frame images for more accurate fiber flux accounting. Requires the defineStripes primitive to be run prior during recipe so INDEX_FIBER and INDEX_ORDER extensions exist to define pixel locations across frame within fiber traces to avoid when finding stray light. Additionally some arm dependent, hard-coded elements are added to the edges of the mask prior to use of photutils.background2D to calculate the background. See [Details](#details-on-primitives-in-primitives_maroonxpy) section of this file for more information.

**End of primitive calls that are run on each stream separately**

- combineFlatStreams: recombines the background-subtracted flat frame data into one processed frame,combining the main (FDDDF illuminated) stream and the 'DFFFD_flats' stream with a simple max comparison at each pix and saves in main stream.

- clearStream: default call to remove the now unneeded non-main stream.

- findStripes and identifyStripes are now called fresh on the fully-illuminated, background-subtracted flat-frame. This is performed for quality assurance as well as correct indexing in meta-data. Outputs are still saved as temporary non-fits-writeable extensions.

- defineFlatStripes: now called on the fully-illuminated, background-subtracted flat-frame with extract_parameter=True. Now saves fiber location info for future science extraction as FITS savable STRIPES_ID and STRIPES_FIBERS. STRIPES_ID and STRIPES_FIBERS contain the by-spectral-order polynomial plate solution for each illuminated fiber that is utilized to define 2D extraction regions in science extractions.

- storeProcessedFlat: standard reduction on the recipe spits out the final frame in ./calibrations/processed_flat/

### SCIENCE FLUX EXTRACTION

- found in default at ./maroonxdr/maroonx/sq/recipes_ECHELLE_SPECT.py
(for makeStripeExtractionCheck info, see test_stripe_retrieval.py test)

This recipe processes MAROON-X echelle spectrum, using previously identified fiber/order traces in a 2D flat, this recipe performs both regular (aka 'box') and optimum extraction to produce 1D extracted spectra for 2D input spectra.
Box extraction is the simple summation of all spatial pixels in a given fiber/order combination. Optimal extraction is per default only applied to fibers illuminated with flat (F) and science (O) input.

Performs checkArm, addDQ, overscanCorrect, and correctImageOrientation identically to their use in the creation of a master dark to ensure all received raw science frames are prepared correctly for use. Given the current implementation of subsequent primitives I think checkArm could be dropped (all frames are independently reduced with individual calls to a reference flat) but I haven't tested this.

- darkSubtraction: Finds the relevant processed_dark frame and creates a DARK_SUBTRACTED extension in the science frame for possible use in the echelle stripe extraction. Currently, the connected processed_dark frame is hardcoded because MX isn't in the caldb.

- extractStripes: creates sparse matrix of each fiber+order given the references previously calculated in the processed flat frame for the BPM, flat, and science frames.  Currently, the connected processed_flat frame is hardcoded because MX isn't in the caldb. Essentially a big wrapper around the _extract_single_stripe function (which is misnamed, is really just retrieval, not extraction) to retrieve all pixels for each element that goes into the optimal extraction.

- optimalExtraction:  automatically performs a normal 2D->1D 'box' extraction and if desired (as is default for target fibers) performs an optimal extraction on each of the sparse echelle spectra found in extractStripes (Requires extractStripes to be run first in recipe). The given corresponding sparse flat field spectra is used to generate normalized 'profiles' that are used as weighting functions in the optimal reduction. The algorithm further checks for cosmic hit outliers and rejects them, iterating the BPM in the process. Output is fits-writeable copies of the (final) BPM, box extraction and error, optimal extraction and error, and a copy of the absolute orders extracted for each fiber. All extensions are always created no matter what is requested to be optimally extracted for ease-of-access following reduction. Essentially a big wrapper around the _optimal_extraction_single_stripe and_box_extract_single_stripe functions to do 2D->1D extraction on all fiber/orders.

- storeProcessedScience: standard reduction on the recipe spits out the final frame in raw science directory.

## DETAILS ON PRIMITIVES IN PRIMITIVES_MAROONX.PY

The following list contains all the primitives in the primitives_maroon.py file in alphabetical order.  These primitives deal with the frames as a whole, as opposed to primitives in primitives_maroonx_echelle.py file, which deal with echelles.

- **addDQ** : This primitive is used to add a DQ extension to the input AstroData object. The value of a pixel in the DQ extension will be the sum of the following: (0=good, 1=bad pixel (found in bad pixel mask), 2=pixel is in the non-linear regime, 4=pixel is saturated). This primitive will trim the BPM to match the input AstroData object(s).
  - Parameters:
    - adinputs : list of AstroData objects with no DQ extensionsuffix: str
            suffix to be added to output files
    - static_bpm: str
            Name of bad pixel mask ("default" -> use default from look-up table)
            If set to None, no static_bpm will be added.
    - user_bpm: str
            Name of the bad pixel mask created by the user from flats and
            darks.  It is an optional BPM that can be added to the static one.
    - illum_mask: bool
            add illumination mask?
  - Returns:
    - adinputs : list of AstroData objects with a DQ extension added to them

- **checkArm** : Check that MX frame arm is consistent through all input files, i.e. BLUE or RED based on data tags. The first file sets the expected value.  Currently, assumes 1 astrodata object comes from 1 single-extension FITS. Need to update if/when original FITS are MEF.
  - Parameters
    - adinputs : list of un-checked MX frames
  - Returns
    - adoutputs : list of frames that have been checked,  always at least the first frame.

- **checkImageOrientation** : Correct image orientation to proper echelle format for MAROON-X.  Flips SCI, if needed, so that left lower corner is bluest wavelength, upper right corner is reddest wavelength.  Resulting echelle orders go from left to right.  MX blue frames start with incorrect orientation for reduction.  This primitive must be called after DQ is established and before any image arithmetic is performed.
  - Parameters
    - adinputs - list of un-checked MX objects
  - Returns
    - adoutputs - same list as inputs, with correct orientation to SCI

- **checkND** : Check that the ND filter on the sim cal fiber is consistent through all MX-input files i.e. illumination is similar intensity as needed for good removal. The first file sets the expected value.
  - Parameters
    - adinputs : list of MX-objects
  - Returns
    - adoutputs : adinputs that pass test

- **addVAR** : Calculates the variance based on the read noise for the chip and the poisson noise(the variance in this case is just the number of photons for each pixel).  The variance is then stored as a FITS extension for each file.
  - Parameters
    - adinputs - list of MX objects without variance extensions
    - read_noise - boolean, whether to include read noise in variance calculations
    - poisson_noise - boolean, whether to include poisson noise in variance calculations
  - Returns
    - adoutputs - list of MX objects with variance extensions
    
- **combineFlatStreams** : Recombines the flat data into one processed frame, combining the main stream pre-processed and the 'source' stream pre-processed with a simple max comparison at each pixel.  Saves in main stream.
  - Parameters
    - 'DFFFD_flats' stream : single MX astrodata object
    - 'main' stream : single MX astrodata object
    - **params needed for access to stream
  - Returns
    - adoutput : single FFFFF_flat MX astrodata object with primary extension data as combined all fiber illuminated flat

- **defineFlatStripes** : Saves fiber location info based on flat field info for stray light removal (extract=False) and for future science extraction (extract=True). Requires the findStripes and identifyStripes primitives to be run prior during recipe so necessary information exists in input extensions. <br> Will remove previous (improperly formatted, but fast) STRIPES_ID and STRIPES_LOC extensions and replace with INDEX_FIBER and INDEX_ORDER pixel map extensions, as needed in straylight removal, and (if extract=True) a FITS savable STRIPES_ID and STRIPES_FIBERS. <br> For a given slit_height, this function extracts the flat field stripes, calculates a box extracted spectrum and normalizes the flat field to generate a 2D pixel map that is used in the straylight removal.
STRIPES_ID and STRIPES_FIBERS contain the by-spectral-order polynomial plate solution for each illuminated fiber that is utilized to define 2D extraction regions in science extractions.
  - Parameters
    - adinputs : single MX astrodata object, is either a DFFFD, FDDDF flat, or combined FFFFF flat
    - slit_height : half pixel height of box in each dimension to
            perform box extraction with
    - extract :  if True, will write STRIPES_ID in fits-acceptable format. Utilized in combined, all fiber illuminated FFFFF_flat

  - Returns
    - adoutput : single MX astrodata object with INDEX_FIBER, INDEX_ORDER extensions and possibly STRIPES_ID and STRIPES_FIBERS extensions

- **findStripes** : Locates and fits stripes in a flat field spectrum. Starting in the central column, the algorithm identifies peaks and traces each stripe to the edge of the detector by following the brightest pixels along each order. It then fits a polynomial to each stripe. To improve algorithm stability, the image is first median filtered and then smoothed with a gaussian. It not only eliminates noise, but also ensures that the cross-section profile of the flat becomes peaked in the middle, which helps to identify the center of each stripe. Choose gauss_filter accordingly. To avoid false positives, only peaks above a certain (relative) intensity threshold are used.
  - Parameters
    - adinputs : single MX astrodata object, is either a DFFFD flat,
            FDDDF flat, or combined FFFFF flat
    - deg_polynomial : degree of the polynomial fit
    - med_filter : median filter parameter
    - gauss_filter_sigma : sigma of the gaussian filter used to
            smooth the image.
    - min_peak : minimum relative peak height
  - Returns
    - adoutput : single MX astrodata object with STRIPES_LOC extension.
        This extension temporarily holds the fits-unsavable fiber information
        before it is utilized and then removed.

- **identifyStripes** : Identifies the stripes by assigning their proper order and fiber number, including correction for the possibilitiy that the spectra have shifted up/down in the cross-dispersion direction since the reference was made.  Requires the findStripes primitive to be run prior during recipe so the stripes are located in the input, i.e. STRIPES_LOC extension exists.
  - Parameters
    - adinputs : single MX astrodata object, is either a DFFFD flat, FDDDF flat, or combined FFFFF flat with STRIPES_LOC extension
    - positions_dir : lookup fits location of nominal y positions and fiber/order labels. Shape is Nx3, columns are [fibers, orders, y_positions], nominally found in lookups/SID
    - selected_fibers : fibers illuminated in the flat, if None, assumes all can work if not given on partially illuminated frame, but best practice is to explicitly identify on function call.
  - Returns
    - adoutput : single MX astrodata object with STRIPES_ID extension.  This extension temporarily holds the fits-unsavable fiber information before it is utilized and then removed. A new extension REMOVED_STRIPES also saves the polynomial info for every stripe that is not identified from the original set inherited from findStripes.

- **removeStrayLight** : Removes stray light from full frame images for more accurate fiber flux accounting. Requires the defineStripes primitive to be run prior during recipe so INDEX_FIBER and INDEX_ORDER extensions exist to define pixel locations across frame within fiber traces to avoid when finding stray light.
  - Parameters
    - adinputs : single MX astrodata object, is either a DFFFD or FDDDF flat
            that has not previously had its stray light removed
    - box_size : pixel height and width of 'mesh_element' used in
            background identification sub-routine
    - filter_size : pixel height and width of window to perform
            background identification sub-routine
    - snapshot : Bool to save difference frame of removed stray light as
            extension STRAYLIGHT_DIFFERENCE
  - Returns
    - adoutput : single MX astrodata object with stray light removed from
            SCI

- **stackDarks** : MX-specific version of stack darks allowing scaling for etalon intensity drift that is in MX 'darks'.  This is needed because although we use them as 'darks' these frames have flux that needs to be stacked appropriately in them. Also utilizes the stackFramesMXCal primitive that scales based on normalizing from average full-frame mean as opposed to the normal stackFrames implementation of first frame's full frame mean (the difference is important for maroonx). Differences to standard implementations are noted with !! comments throughout.
  - Parameters are inherited from stackFramesMXCal
  - Returns
    - adoutputs : Master Dark. This list contains only one element. The list format is maintained so this primitive is consistent with all the others.

- **stackFlats** : stacks based on the stackFramesMXCal primitive that scales based on average full-frame mean.
  - Parameters are inherited from stackFramesMXCal
  - Returns
    - adoutputs : Master Flat. This list contains only one element. The list format is maintained so this primitive is consistent with all the others.

- **stackFramesMXCal** : MX-specific version of stackFrames for calibration frames - changes scaling to average full frame mean to purposely scale by etalon flux and its drift between calibration exposures, this function should not be used to combine MX science frames.
  - Parameters
    - adinputs : List of AstroData objects to be combined.
    - suffix : str
            Suffix to be added to output files.
    - apply_dq : bool
            Apply DQ mask to data before combining?
    - nlow, nhigh : int
            Number of low and high pixels to reject, for the 'minmax' method.
            The way it works is inherited from IRAF: the fraction is specified
            as the number of  high  and low  pixels,  the  nhigh and nlow
            parameters, when data from all the input images are used.  If
            pixels  have  been  rejected  by offseting,  masking, or
            thresholding then a matching fraction of the remaining pixels,
            truncated to an integer, are used.  Thus::

                nl = n * nlow/nimages + 0.001
                nh = n * nhigh/nimages + 0.001

            where n is the number of pixels  surviving  offseting,  masking,
            and  thresholding,  nimages  is the number of input images, nlow
            and nhigh are task parameters  and  nl  and  nh  are  the  final
            number  of  low  and high pixels rejected by the algorithm.  The
            factor of 0.001 is to adjust for rounding of the ratio.
    - operation : str
            Combine method.
    - reject_method : str
            Pixel rejection method (none, minmax, sigclip, varclip).
    - zero : bool
            Apply zero-level offset to match background levels?
    - scale : bool
            Scale images to the same intensity?
    - memory : float or None
            Available memory (in GB) for stacking calculations.
    - statsec : str
            Section for statistics.
    - separate_ext : bool
            Handle extensions separately?

  - Returns
    - adoutputs : Sky stacked image. This list contains only one element. The list format is maintained so this primitive is consistent with all the others.

  - Raises
    - IOError
            If the number of extensions in any of the `AstroData` objects is
            different.
    - IOError
            If the shape of any extension in any `AstroData` object is different.
    - AssertError
            If any of the `.gain()` descriptors is None.
    - AssertError
            If any of the `.read_noise()` descriptors is None.

- **validateData** : MAROON-X-specific version of validateData to ignore the invalid WCS exception.
  - Parameters
    - adinputs : List of unchecked AstroData objects.
  - Returns
    - adinputs : List of checked AstroData objects.

## DETAILS ON PRIMITIVES IN PRIMITIVES_MAROONX_ECHELLE.PY

This section contains details on the primitives in primitives_maroonx_echelle.py. The primitives are listed in alphabetical order.

- _box_extract_single_stripe: Box extraction of a single stripe.
  - Parameters
    - stripe (sparse.matrix): stripe to be extracted
    - mask (np.ndarray): binary pixel mask. Where mask==1: values will contribute. Same shape as stripe.toarray()
  - Returns:
    - np.ndarray: box extracted flux

- darkSubtraction: Finds the dark frame in association with the adinput and creates a dark subtracted extension that can be requested during stripe extraction
  - Parameters
    - adinputs: AstroData object(s) for which dark subtraction is to be performed
    - dark: (optional) adinput of relevant processed dark
    - individual: (bool) if True, creates a calib call for each individual science frame.
            If False, groups frames into exposure time and ND_filter and calls one calib per group
            based on first frame in group.
  - Returns
    - adinputs with additional image extension of dark subtracted full frame.

- _extract_single_stripe:  Extracts single stripe from 2d image.  This function returns a sparse matrix containing all relevant pixels for a single stripe for a given polynomial p and a given slit height.  This is a helper function for extractStripes, and should be altered with caution.
  - Parameters
    - polynomials (np.ndarray): polynomial coefficients
    - img (np.ndarray): 2d echelle spectrum
    - slit_height (int): total slit height in pixel to be extracted
  - Returns
    - scipy.sparse.csc_matrix: extracted spectrum

- extractStripes: Extracts the stripes from the original 2D spectrum to a sparse array, containing only relevant pixels. <br> This function marks all relevant pixels for extraction.  Reinterpreting the flat reference it iterates over all stripes in the image and saves a sparse matrix for each stripe.
  - Parameters
    - adinputs
    - flat: adinput of relevant processed flat, as processed, will have the STRIPES_ID and STRIPES_FIBERS extensions needed
    - dark: (optional) adinput of relevant processed dark
    - skip_dark: if dark given, which individual fibers dark subtraction should be skipped
    - slit_height (int): total slit height in px
    - test_extraction (bool): used in unit test for this function, saves science extraction, flat extraction, and the bpm-extraction in FITS-readable format (STRIPES, F_STRIPES, STRIPES_MASK)
    - individual: (bool) if False uses one calib call for all frames per arm, if True performs a calib call for each frame
  - Returns
    - adinputs with sparse matrices added holding the 2D extractions for each fiber/order for the science frame, flat frame, and BPM (STRIPES, F_STRIPES, STRIPES_MASK)if test_extraction==True, the extractions are FITS-readable and not sparse matrix format

- optimalExtraction: Optimal extraction of the 2d echelle spectrum.<br>  This function performs an optimal extraction of a 2d echelle spectrum. A given flat field spectrum is used to generate normalized 'profiles' that are used as weighting functions for the spectrum that is going to be extracted.  The algorithm further checks for outliers and rejects them.  This is to prevent contributions from cosmic hits.
  - Parameters
    - adinputs with STRIPES, F_STRIPES, and STRIPES_MASKS 'extensions' as dicts of sparse arrays
    - opt_extraction (list): fibers considered for optimal extraction
    - back_var (float): manual background variance for frame
    - full_output (bool): if True, an additional set of intermediate products will be returned / saved penalty (float): scaling penalty factor for mismatch correction between flat field profile and science spectrum during optimal extraction
    - s_clip (float): sigma-clipping paramter during optimal extraction
  - Returns
    - adinputs with optimal and box extracted orders for each fiber as well as uncertainties and the bad pixel mask result from the optimal extraction

- _optimal_extraction_single_stripe: Performs optimal extraction of a single stripe.  Based on the algorithm described by Horne et al. 1986, PASP, 98, 609.  This is a helper function for optimalExtraction and should be altered with caution.
  - Parameters
    - stripe (scipy.sparse.spmatrix): sparse matrix stripe from science
                frame to be extracted
    - flat_stripe (scipy.sparse.spmatrix): sparse matrix stripe from flat frame to be used as profile
    - gain (float): detector gain factor (conversion photons -> DN, given in e-/DN )
    - read_noise (float): typical detector read noise given as the variance in DN
    - back_var (np.ndarray): background variance (from scattered light model or dark, otherwise 0)
    - mask (np.ndarray): bad pixel mask. Note: the mask will be modified by iterative algorithm that looks for outliers during optimal extraction
    - full_output (bool): if True, returns all intermediate results for debugging/testing
    - s_clip (float): sigma clipping value for optimal extraction penalty (float): scaling factor for global per-order profile mismatch correction between expected (flat) and found (science).  Set 0 for no correction
  - Returns
    - tuple(np.ndarray, np.ndarray, dict): (optimal extracted spectrum, box extracted spectrum, dict of additional intermediate results if full_output was True)

## TESTS

The tests are divided into three categories: echelle tests, image tests, and complete tests.  The echelle tests are unit tests for the primitives in primitives_maroonx_echelle.py, while image tests are unit tests for the primitives in primitives_maroonx.py. These tests are implemented using pytest, so please make sure you have pytest installed in your conda environment. Complete tests on the other hand test the entire reduction process for darks, flats and science files, and are not implemented using pytest.  The end user should only have to worry about the complete tests, and should not have to worry about the echelle or image tests. 
### COMPLETE TESTS

found in ./maroonxdr/maroonx/tests/complete/

The complete tests are not written using pytest.  Instead these are just 3 simple scripts that allow the end user to test the entire reduction process for darks, flats, and science files.  The end user should only have to worry about these tests, and should not have to worry about the echelle or image tests.  The suser should be able to run these tests by simply running python xxxx_test.py in the terminal.
### IMAGE TESTS

found in ./maroonxdr/maroonx/tests/image/

To ensure that the unit tests work, we suggest putting your fits files in a folder called science_dir inside the MAROONXDR directory. If you put your fits files in a different directory, you will need to change the path in the test files.

#### KNOWN ISSUES

Currently, 2 tests in checkND and the tests in findStripes fail.  This is expected and will be investigated later.

#### TESTS IN ALPHABETICAL ORDER

- test_file_sorting.py: Composed of three tests functions that test the checkArm, separateFlatStream, and combineFlatStream primitives respectively for use cases.
  - test_checkArm_collection_and_rejection: is purposely given frames of both the red and blue arms to ensure that the primitive tested correctly warns and truncates output to just the type that the first given frame has.
    - Parameters
      - caplog : fixture
      - filename_r : str
      - filename_b : str
    - Returns
      - None

  - test_separating_flat_streams: Given a series including both types of raw illuminated flat frames of a single arm, the primitive is tested on whether it correctly separates the sets and creates the second stream against the directly calculated expectation in the test.
    - Parameters
      - caplog : fixture
      - filename_r : str
      - filename_b : str
    - Returns
      - None

  - test_combining_flat_streams: Given two stacked flat frames of a given arm and the resulting 'combined' frame that the primitive creates and tests this directly against the directly calculated expectation in the test.
    - Parameters
      - caplog : fixture
      - filename_r : str
      - filename_b : str
    - Returns
      - None

- test_image_orientation_corrector.py: Composed of two test functions that tests the correctImageOrientation primitive for its reaction to red frames and to blue frames.
  - test_correctImageOrientation_does_not_change_red_frames: Given a red frame, the test function checks that the primitive output pixel values remain the same as they started (i.e. no transformation is made on the data)
    - Parameters
      - caplog : fixture
      - filename : str
    - Returns
        - None

  - test_correctImageOrientation_flips_blue_frames: Given a blue frame, the test function checks that the primitive output pixel values did undergo flips in both the vertical and horizontal directions.
    - Parameters
      - caplog : fixture
      - filename : str
    - Returns
      - None

- test_ND_filter_check.py: Three test functions that check the rejection capabilities of the checkND primitive.
  - test_nd_filter_good_series: Given a series of frames that contain the same ND filter value, the primitive output is checked to see that it holds them
    - Parameters
      - caplog : fixture
      - filename : str
    - Returns
      - None
  - test_nd_filter_subgood_series: Given a series of frames that hold some frames of the same ND filter value (equal to the first frame's), the primitive output is checked to see that the undesired frames no longer remain, the frames of the same ND filter do remain, and that a warning was given.
    - Parameters
      - caplog : fixture
      - filename : str
    - Returns
      - None
  - test_nd_filter_bad_series: given a series of frames of which the first frame has a unique ND filter value, the primitive is checked for flagging the specific 'first frame uniqueness' IO error
    - Parameters
      - caplog : fixture
      - filename : str
    - Returns
      - None

- test_stray_light_removal.py: a single test file that tests the data manipulation of the stray light removal element of the flat creation recipe against a file previously stored using the makeStrayLightCheck non-default flat reduction.<br> By default the stray light removal step values are not stored while performing flat creation - the makeStrayLightCheck uses the snapshot=True parameter in the removeStrayLight primitive to save a fits-writeable extension of the stray light calculation in each partially illuminated flat and the reduction is stopped and stored at that point. <br> Afterwords that previously calculated data is used to test the continued success of the removeStrayLight primitive by sending the frames, that still contain their raw data, through the steps of the normal recipe and then checking that the newly calculated difference is the same as the previously stored difference in the unique extension that remains unaffected by the new partial reduction.
  - Parameters
    - caplog : fixture
    - filename_r : str
    - filename_b : str
  - Returns
    - None

- test_stripe_finding.py: three tests and a private method that directly test the continued success of the findStripes, identifyStripes, and defineFlatStripes methods.
  - test_full_stripe_definition: given a previously reduced masterflat frame run findStripes, identifyStripes, and defineFlatStripes on the 2D image and assert that for every fiber/order found in the original reduction it is found again by the primitive.
    - Parameters
      - caplog : fixture
      - filename : str

    - Returns
      - None
  - test_identify_stripes: given a previously reduced masterflat frame run findStripes and  identifyStripes on the 2D image and test that for every fiber/order previously found, its trace is identified again as the same real fiber/oder by the identifyStripes primitive. Also that for any other traces newly identified by the findStripes primitive that identifyStripes correctly has found them to be not real fiber/orders as defined by their exclusion from those previously saved as well as their id values following new identification.
    - Parameters
      - caplog : fixture
      - filename : str

    - Returns
      - None
  - test_find_stripes: given a previously reduced masterflat frame run findStripes. Test that the traces that newly found include exactly all previously saved fiber/order traces combined with those previously rejected during the original reduction.
    - Parameters
      - caplog : fixture
      - filename : str

    - Returns
      - None

### ECHELLE TESTS

found in ./maroonxdr/maroonx/tests/echelle_extraction/

- test_dark_subtraction.py: Composed of only one test, 