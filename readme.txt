INTRODUCTION
The red-optical, fiber-fed echelle-spectrograph has no movable parts and only one mode. At each observed wavelength the 5 fibers are arranged along the cross-dispersion direction. The middle three fibers are pupil-sliced fractions of the on-sky target fiber, the outer two fiber traces are an off-target on-sky fiber and a fiber known as the 'sim. cal. fiber' which can be fed calibration light during science frames (as well as during other frames). All observed light is split with a dichroic to a ‘blue arm’ (491-670nm) and a ‘red arm’ (649-920nm). Both arms terminate in CCD detectors with 16-bit 4400x4400 arrays, but the 'blue' detector has 4 reads and the 'red' detector has 2. 

The repo is currently a re-implementation of the 'flux extraction' part of the internal U-Chicago pipeline. This first step occurs before a second wavelength calibration step of the internal two-step pipeline.

Broadly the repo currently makes a few products that are recognizable as standard products:
1) master darks
2) master flats
3) 2D echelle optimally extracted reduced science frames


==============================
RAW FRAME EXPECTATIONS
Repo currently expects on any raw MAROON-X fits read:

1) all raw maroon-x fits files to have keyword INSTRUME = 'MAROON-X' (using the matches_data method repo then internally removes hyphen in instrument definition to agree with repo names)

2) all raw maroon-x fits files with corresponding red and blue arm exposures as completely separate files (currently based on file name in data_tag: _tag_arm)

3) All raw fits frames are full-frame images with no changing bias/data/array/detector section values. Each read section is treated as a separate data section with a unique name. The values are stored in the ./maroonx_instruments/maroonx/lookups.py file.

4) all raw maroon-x fits files to have defined HIERARCH FIBER1 through FIBER5 for all frames (multiple data_tags utilize these entries to to set the associated TagSet for any file type)
-- currently utilzed data_tags are _tag_dark, _tag_flat, and _tag_science which currently utilize 'Dark'+'Etalon', 'Flat', and 'OBJECT'+'ETALON' entry values respectively to make 'DARK'+'CAL', 'FLAT'+'CAL', and 'SCI' TagSet associations respectively.

5) Currently with any raw frame, HIERARCH MAROONX IMAGE ORIENTATION HORIZONTAL FLIP and HIERARCH MAROONX IMAGE ORIENTATION VERTICAL FLIP are utilized to re-orient blue frames to have the same dispersion and cross-dispersion that red frames have. This is redundant with _tag_arm.
 
=============================
RUN OPTIONS
Current run options:
Performing a standard 'reduce' command on all frames of a particular arm and file type can be achived through the standard methods

e.g.
all_files = glob.glob('/some_path/*_r_*.fits')
all_files.sort()
all_flats = dataselect.select_data(all_files, ['FLAT'])
myreduce = Reduce()
myreduce.files.extend(all_flats)
myreduce.drpkg = 'maroonxdr'
myreduce.runr()

=============================
MASTER DARK CREATION DETAILS
- found at ./maroonxdr/maroonx/sq/recipes_DARK.py
This recipe performs the standardization and corrections needed to convert the raw input dark images into a single processed dark image. This output processed dark is stored on disk using storeProcessedDark and has a name equal to the name of the first input bias image with "_dark.fits" appended. The background in an un-illuminated frame is very low for exposure times of less than 900s and likely doesn't warrant a dark subtraction, however most science frames are taken with the simultaneous calibration fiber illuminated with the FP etalon. The extended wings of the etalon reach into one of the science fibers with a few 10 counts. To remove these wings and any broad diffuse (illumination independent) background, DDDDE frames are taken in daytime to construct a DDDDE processed dark. These darks are specific for different exposure times (i.e. ND filter settings) and should be taken close in time (within a day or two) to the science frame as the etalon source brightness can be time variable.

notes on recipe in chronological order of primitive calls:

checkArm: Robust to being fed red and blue arm flats with the checkArm primitive:
will only use the flats in agreement with the first dark's arm and will log warning, error if only one

checkND: Robust to being fed flats at different exposure times with the checkND primitive:
will only use the flats in agreement with the first dark's ND and will log warning, error if only one

addDQ: Adds BPM in ./maroonxdr/maroonx/lookups/BPM/ depending on data sereis arm with the addDQ overwrite of standard call in maroonx primitives (only needed because we are not in caldb). Differences are noted with !! commments

overscanCorrect: Overscan corrects each frames depending on arm and instrument lookup values.

correctImageOrientation: Corrects image orientation, as needed with blue arm, so that left lower corner is bluest wavelength, upper right corner is reddest wavelength. As implemented, this goes (as needed) after overscan correction because it changes blue arm pixel index identification. This consistency of orientation is needed for steps in the science reduction and wavelength association.

stackDarks: Stacks darks using maroonx primitive overwrite of standard. This is needed because although we use them as 'darks' these frames have flux that needs to be stacked appropriately in them. Also utilizes the stackFramesMXCal primitive that scales based on normalizing from average full-frame mean as opposed to the normal stackFrames implementation of first frame's full frame mean (the difference is important for maroonx). Differences to standard implementations are noted with !! comments throughout.

storeProcessedDark: Standard reduction on the recipe spits out the final frame in ./calibrations/processed_dark/


=============================
MASTER FLAT CREATION DETAILS
- found in default at ./maroonxdr/maroonx/sq/recipes_FLAT_SPECT.py
(for makeStrayLightCheck info, see test_stray_light_removal.py test)
This recipe performs the standardization and corrections needed to convert the raw input flat images into a single stacked flat image. This output processed flat is stored on disk using storeProcessedFlat and has a name equal to the name of the first input bias image with "_FFFFF_flat.fits" appended. A processed flatfield is required to perform optimal flux extraction and to determine the blaze function. Due to the stability of the spectrograph, one processed flatfield is typically 'valid' for at least a two-week period. Cross-comparison between different processed flatfields taken months apart have not been conducted so far. The recipe expects two sets of sub-illumated flat frames (as needed for best background subtraction processing).

notes on recipe in chronological order of primitive calls:

Performs checkArm, checkND, addDQ, overscanCorrect, and correctImageOrientation identically to their use in the creation of a master dark to ensure all recieved raw flat frames have identical setup and are prepared correctly for use. 

seperateFlatStreams: uses fiber_setup data_descriptor to identify which of the given frames have flat illumation in the non-target fibers (known internally in the code as FDDDF) and stay in the main stream and which have flat illumination in the target fibers (i.e. DFFFD) that get moved to the newly created DFFFD_flats stream. Warns if no frames are found in either of the two patterns.

---starting here primitive calls are run on both the main stream and the 'DFFFD_flats' stream

stackFlats: stacks based on the stackFramesMXCal primitive that scales based on average full-frame mean.

findStripes: locates and fits stripes in the stacked flat field spectrum. Starting in the central column, the algorithm identifies peaks and traces each stripe to the edge of the detector by following the brightest pixels along each order. It then fits a polynomial to each stripe. This information is saved in a non-fits-writeable extension of the astrodata object. I.E. find all fiber traces in the 2D frame, save the polynomal that defines their pixel locations.

identifyStripes: identifies the stripes to their proper order and fiber number, including correction for the possibilitiy that the spectra have shifted up/down in the cross-dispersion direction since the reference was made. Reference fits tables are found in ./maroonxdr/maroonx/lookups/SID/ using the internal _get_sid_filename primitive function. This primitive requires the findStripes primitive to be run prior during recipe so the stripes pixel traces are located in the input, i.e. STRIPES_LOC extension exists. The proper order and fiber number identification is also held in a non-fits-writeable extension. I.E. decide which traces are real fiber/orders and zip the pixel location polynomials with the fiber/order IDs.

defineFlatStripes: The first two times this primitive is called (extract==Default==False) is to map all pixels in each stacked flat field frame that have been identified as part of a real fiber trace so that they can be masked in the following background modelling. Requires the findStripes and identifyStripes primitives to be run prior during recipe so necessary information exists in input extensions. Will remove previous (improperly formatted, but fast) STRIPES_ID and STRIPES_LOC extensions and replace with INDEX_FIBER and INDEX_ORDER pixel map extensions, as needed in straylight removal. Previous, temporary non-fits-writeable extensions are dropped. I.E. with extract == False creates two 2D composite map frames with the pixel entries being their fibers and orders respectively for all pixels found (by identifyStripes) to be in a fiber trace.

removeStrayLight: Removes stray light from full frame images for more accurate fiber flux accounting. Requires the defineStripes primitive to be run prior during recipe so INDEX_FIBER and INDEX_ORDER extensions exist to define pixel locations across frame within fiber traces to avoid when finding stray light. Additionally some arm dependent, hard-coded elements are added to the edges of the mask prior to use of photutils.background2D to calculate the background. 

---end of primitive calls that are run on eacch stream separately

combineFlatStreams: recombines the background-subtracted flat frame data into one processed frame,combining the main (FDDDF illumated) stream and the 'DFFFD_flats' stream with a simple max comparison at each pix and saves in main stream.

clearStream: default call to remove the now unneeded non-main stream.

findStripes and identifyStripes are now called fresh on the fully-illumated, background-subtracted flat-frame. This is performed for quality assurance as well as correct indexing in meta-data. Outputs are still saved as temporary non-fits-writeable extensions.

defineFlatStripes: now called on the fully-illumated, background-subtracted flat-frame with extract_parameter=True. Now saves fiber location info for future science extraction as FITS savable STRIPES_ID and STRIPES_FIBERS. STRIPES_ID and STRIPES_FIBERS contain the by-spectral-order polynomial plate solution for each illuminated fiber that is utilized to define 2D extraction regions in science extractions.

storeProcessedFlat: standard reduction on the recipe spits out the final frame in ./calibrations/processed_flat/
        
    
=============================
SCIENCE FLUX EXTRACTION DETAILS    
- found in default at ./maroonxdr/maroonx/sq/recipes_ECHELLE_SPECT.py
(for makeStripeExtractionCheck info, see test_stripe_retrieval.py test)
This recipe processes MAROON-X echelle spectrum, using previously identified fiber/order traces in a 2D flat, this recipe performs both regular (aka 'box') and optimum extraction to produce 1D extracted spectra for 2D input spectra.
Box extraction is the simple summation of all spatial pixels in a given fiber/order combination. Optimal extraction is per default only applied to fibers illuminated with flat (F) and science (O) input.        
        
Performs checkArm, addDQ, overscanCorrect, and correctImageOrientation identically to their use in the creation of a master dark to ensure all recieved raw science frames are prepared correctly for use. Given the current implementation of subsequent primitives I think checkArm could be dropped (all frames are indepenedently reduced with individual calls to a reference flat) but I haven't tested this.

darkSubtraction: Finds the relevant processed_dark frame and creates a DARK_SUBTRACTED extension in the science frame
for possible use in the echelle stripe extraction. Currently, the connected processed_dark frame is hardcoded because MX isn't in the caldb.

extractStripes: creates sparse matrix of each fiber+order given the references previously calculated in the processed flat frame for the BPM, flat, and science frames.
Currently, the connected processed_flat frame is hardcoded because MX isn't in the caldb. Essentially a big wrapper around the _extract_single_stripe function (which is misnamed, is really just retrieval, not extraction) to retrieve all pixels for each element that goes into the optimal extraction.

optimalExtraction:  automatically performs a normal 2D->1D 'box' extraction and if desired (as is default for target fibers) performs an optimal extraction on each of the sparse echelle spectra found in extractStripes (Requires extractStripes to be run first in recipe). The given corresponding sparse flat field spectra is used to generate normalized 'profiles' that are used as weighting functions in the optimal reduction. The algorithm further checks for cosmic hit outliers and rejects them, iterating the BPM in the process. Output is fits-writeable copies of the (final) BPM, box extraction and error, optimal extraction and error, and a copy of the absolute orders extracted for each fiber. All extensions are always created no matter what is requested to be optimally extracted for ease-of-access following reduction. Essentially a big wrapper around the _optimal_extraction_single_stripe and _box_extract_single_stripe functions to do 2D->1D extraction on all fiber/orders.
        
storeProcessedScience: standard reduction on the recipe spits out the final frame in raw science directory. 

=============================
IMAGE TESTS
found in ./maroonxdr/maroonx/tests/image/

test_file_sorting.py: composed of three tests functions that test the checkArm, seperateFlatStream, and combineFlatStream primitives respectively for use cases.
-test_checkArm_collection_and_rejection: is purposely given frames of both the red and blue arms to ensure that the primitive tested correctly warns and truncates output to just the type that the first given frame has
-test_separating_flat_streams: given a series including both types of raw illuminated flat frames of a single arm, the primitive is tested on whether it correctly separates the sets and creates the second stream against the directly calculated expectation in the test
-test_combining_flat_streams: is given two stacked flat frames of a given arm and the resulting 'combined' frame that the primitive creates is tested directly against the directly calculated expectation in the test.


test_image_orientation_corrector.py: composed of two test functions that tests the correctImageOrientation primitive for its reaction to red frames and to blue frames.
-test_correctImageOrientation_does_not_change_red_frames: given a red frame, the test function checks that the primitive output pixel values remain the same as they started (i.e. no transformation is made on the data)
-test_correctImageOrientation_flips_blue_frames: given a blue frame, the test function checks that the primitive output pixel values did undergo flips in both the vertical and horizontal directions.


test_ND_filter_check.py: three test functions that check the rejection capabilites of the checkND primitive.
-test_nd_filter_good_series: given a series of frames that contain the same ND filter value, the primitive output is checked to see that it holds them
-test_nd_filter_subgood_series: given a series of frames that hold some frames of the same ND filter value (equal to the first frame's), the primitive output is checked to see that the undesired frames no longer remain, the frames of the same ND filter do remain, and that a warning was given.
-test_nd_filter_bad_series: given a series of frames of which the first frame has a unique ND filter value, the primitive is checked for flagging the specific 'first frame uniqueness' IO error 


test_stray_light_removal.py: a single test file that tests the data manipulation of the stray light removal element of the flat creation recipe against a file previously stored using the makeStrayLightCheck non-default flat reduction. 
By default the stray light removal step values are not stored while performing flat creation - the makeStrayLightCheck uses the snapshot=True parameter in the removeStrayLight primitive to save a fits-writeable extension of the stray light calculation in each partially illuminated flat and the reduction is stopped and stored at that point. 
Afterwords that previously calculated data is used to test the continued success of the removeStrayLight primitive by sending the frames, that still contain their raw data, through the steps of the normal recipe and then checking that the newly calculated difference is the same as the previously stored difference in the unique extension that remains unaffected by the new partial reduction.


test_stripe_finding.py: three tests and a private method that directly test the continued success of the findStripes, identifyStripes, and defineFlatStripes methods. 

-test_full_stripe_definition: given a previously reduced masterflat frame run findStripes, identifyStripes, and defineFlatStripes on the 2D image and assert that for every fiber/order found in the original reduction it is found again by the primitive.

-test_identify_stripes: given a previously reduced masterflat frame run findStripes and  identifyStripes on the 2D image and test that for every fiber/order previously found, its trace is identified again as the same real fiber/oder by the identifyStripes primitive. Also that for any other traces newly identified by the findStripes primitive that identifyStripes correctly has found them to be not real fiber/orders as defined by their exclusion from those previously saved as well as their id values following new identification. 

-test_find_stripes: given a previously reduced masterflat frame run findStripes. Test that the traces that newly found include exactly all previously saved fiber/order traces combined with those previously rejected during the original reduction.


=============================
ECHELLE EXTRACTION TESTS
found in ./maroonxdr/maroonx/tests/echelle_extraction/
*issue with standard reduction_flat and reduction_dark values being to short to save full filename for proper use to make sure test primitive uses same file as previous reduction*

test_dark_subtraction.py: tests that the re-creation of a dark subtracted science frame with the darkSubtraction primitive is as it was before with comparison to the actual pixel values in the
DARK_SUBTRACTED extension.

test_stripe_retrieval.py: tests the creation of the flat, science frame, and BPM, sparse matricies for all fiber/orders by the extractStripes primitive given a science frame that has been previously partially reduced using the makeStripeExtractionCheck non-default echelle_spect recipe. Tests both red arm and blue arm frame-extraction as currently written. 
The makeStripeExtractionCheck runs all primitives of the regular recipe from a raw input science frame and uses the 'test_extraction' parameter = True to have the sparse matricies all be rewritten in (dense) fits-writeable format extensions.
In the test this previously calculated information is all np.testing.assert_allclose against the freshly calculated values based on the primitive running during the test on the 2D image.


test_extraction.py: For both a red arm frame and a blue arm frame, test the optimal extraction primitive results against previously run results as saved in an extraction-complete science frame.

 
