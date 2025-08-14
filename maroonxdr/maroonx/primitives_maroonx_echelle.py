"""
Primitives for MAROON-X echelle data reduction.  Primities in this file are
focused on the code to produce 1D extracted spectra from 2D spectra
(as compared to primitives_maroonx_2D, which is focused on primitives
that focus on the 2D spectra themselves).
"""
# ------------------------------------------------------------------------------

# from geminidr.gemini.lookups import DQ_definitions as DQ
import copy
from pathlib import Path
import os

import astrodata
import matplotlib.pyplot as plt
import numba
import numpy as np
from geminidr.core import Spect
from gempy.adlibrary import dataselect
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override
from scipy import sparse
from scipy.ndimage import median_filter

from . import parameters_maroonx_echelle
from .primitives_maroonx_2D import MAROONX

# ------------------------------------------------------------------------------

def _get_calibration_flat_path():
    """
    Get the path for the calibration flat file.
    Should probably be deprecated when dragons calib is implemented.
    """
    cwd = Path(os.getcwd())
    return cwd / 'calibrations' / 'processed_flat'

def _get_calibration_dark_path():
    """
    Get the path for the calibration dark file.
    Should probably be deprecated when dragons calib is implemented.
    """
    cwd = Path(os.getcwd())
    return cwd / 'calibrations' / 'processed_dark'

def _get_calibration_dark_coeff(arm_tag):
    """
    Match and return calibration dark coefficient file as astrodata object.
    Should probably be deprecated when dragons calib is implemented.
    """
    # Get the calibration darks
    calib_dark_path = _get_calibration_dark_path()
    files = dataselect.select_data(
        list(calib_dark_path.glob('*.fits')), tags=['DARK_COEFF', arm_tag])

    if len(files) == 0:
        raise ValueError(f'No dark coefficient file found for arm tag: {arm_tag}')
    elif len(files) > 1:
        raise ValueError(f'Multiple calibration darks found for arm tag: {arm_tag}')
    return astrodata.open(files[0])

# def _get_calibration_dark(adinputs):
#     """
#     Match and return calibration darks as astrodata list.
#     Should probably be deprecated when dragons calib is implemented.
#     """
#     # Get the calibration darks
#     calib_dark_path = _get_calibration_dark_path()
#     if "SCI" in adinputs[0].tags:
#         # need to create a synthetic dark
#     elif 


@parameter_override
class MAROONXEchelle(MAROONX, Spect):
    """
    This class contains primitives that applies to all MAROON-X echelle
    data.  Specifically, this class is focused on the code to produce 1D
    extracted spectra from 2D spectra.   Additional steps to produce wavelength
    calibrations are in the MaroonXSpectrum class.
    """

    tagset = {'GEMINI', 'MAROONX', 'ECHELLE'}

    def __init__(self, adinputs, **kwargs):
        super(MAROONXEchelle, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'maroonxdr.maroonx.lookups'
        self._param_update(parameters_maroonx_echelle)

    def _getSyntheticDark(self, adinputs=None, **params):
        """
        Get the path for the synthetic dark file.
        Should probably be deprecated when dragons calib is implemented.
        """
        # adouts = []
        # for ad in adinputs:
        #     # get the coefficients ad
        #     # calculate synthetic dark
        #     # append to adouts
        return

    def attachSyntheticDark(self, adinputs=None, dark_coeff=None, individual=False, **params):
        """
        Finds or creates synthetic dark frames and attaches them to adinputs as SYNTH_DARK extension.
        
        Creates synthetic darks with caching behavior controlled by individual parameter.
        
        Parameters
        ----------
        adinputs: AstroData object(s) to attach darks to
        dark_coeff: (optional) adinput of dark coefficients file
        individual: (bool) if True, creates unique dark for each frame (no reuse).
            If False, reuses darks for frames with same exposure time, ND filter, and arm.
            
        Returns
        -------
        adinputs with SYNTH_DARK extension added to each frame
        """
        log = self.log
        #log.debug(gt.log_message("primitive", self.myself(), "starting"))
        #timestamp_key = self.timestamp_keys[self.myself()]
        
        # Cache for storing created darks
        dark_cache = {}
        
        for ad in adinputs:

            if hasattr(ad[0], 'SYNTH_DARK'):
                log.fullinfo(f"Skipped. SYNTH_DARK extension found in {ad.filename}.")
                continue

            exptime = ad.exposure_time()
            nd_filter = ad.filter_orientation()['ND']
            arm_tag = 'BLUE' if 'BLUE' in ad.tags else 'RED'
            
            # Create cache key based on individual parameter
            if individual:
                cache_key = ad.filename  # Unique key - no reuse
            else:
                cache_key = (exptime, nd_filter, arm_tag)  # Grouping key - allows reuse
            
            if cache_key not in dark_cache:
                # Create new synthetic dark
                if dark_coeff is None:
                    dark_coeff_ad = _get_calibration_dark_coeff(arm_tag)
                else:
                    dark_coeff_ad = dark_coeff
                
                if dark_coeff_ad is not None:
                    synthetic_dark, actual_exptime, actual_nd, factor = create_synthetic_dark(
                        dark_coeff_ad, exptime_value=exptime, nd_value=nd_filter
                    )
                    dark_cache[cache_key] = synthetic_dark
                    log.fullinfo(f"Created synthetic dark for {ad.filename}: "
                                f"exptime={actual_exptime:.1f}s, nd={actual_nd:.2f}, factor={factor:.2f}")
                else:
                    log.warning(f"No dark coefficients found for {arm_tag} arm, {ad.filename}")
                    dark_cache[cache_key] = None
            
            # Attach the cached dark to this frame
            if dark_cache[cache_key] is not None:
                ad[0].SYNTH_DARK = dark_cache[cache_key].copy()
                # gt.mark_history(ad, primname=self.myself(),
                #             keyword='SYNTH_DARK',
                #             comment=f"synthetic_{arm_tag.lower()}_{exptime}s_{nd_filter}")
            else:
                ad[0].SYNTH_DARK = None
                log.warning(f"No synthetic dark attached to {ad.filename}")
        
        # gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def attachDarkSubtraction(self, adinputs=None, **params):
        """
        Performs dark subtraction using the SYNTH_DARK extension to create DARK_SUBTRACTED extension.
        
        This primitive assumes that attachDark() has already been run to create the SYNTH_DARK extension.
        It performs pixel-by-pixel subtraction: DARK_SUBTRACTED = data - SYNTH_DARK
        
        Parameters
        ----------
        adinputs: AstroData object(s) with SYNTH_DARK extension
        
        Returns
        -------
        adinputs with DARK_SUBTRACTED extension added
        """
        log = self.log
        # log.debug(gt.log_message("primitive", self.myself(), "starting"))
        # timestamp_key = self.timestamp_keys[self.myself()]
      
        dark_type = params.get('dark_type')
        
        if dark_type == 'synthetic':
            self.attachSyntheticDark()
        else:
            # self.attachClosestDark()
            raise NotImplementedError("Only synthetic dark type is supported")

        adoutputs = []
        
        for ad in adinputs:
            adout = copy.deepcopy(ad)
            
            
            # Check if SYNTH_DARK extension exists
            if hasattr(ad[0], 'SYNTH_DARK') and ad[0].SYNTH_DARK is not None:
                # Perform pixel-by-pixel subtraction
                adout[0].DARK_SUBTRACTED = ad[0].data - ad[0].SYNTH_DARK
                log.fullinfo(f"Dark subtraction completed for {ad.filename}")
                
                # Mark history
                # gt.mark_history(adout, primname=self.myself(),
                #             keyword='DARK_SUBTRACTED',
                #             comment="synthetic_dark_subtracted")
            else:
                log.warning(f"Dark subtraction skipped. "
                        f"No SYNTH_DARK extension found for {ad.filename}.")
                # Copy data without subtraction
                # adout[0].DARK_SUBTRACTED = ad[0].data.copy()
            
            adoutputs.append(adout)
        
        # gt.mark_history(adoutputs, primname=self.myself(), keyword=timestamp_key)
        return adoutputs

    def darkSubtraction_old(self, adinputs=None, dark=None, individual=False,
                        **params):
        """
        Finds the dark frame in association with the adinput and creates a
        dark subtracted extension that can be requested during stripe extraction

        TODO: In current format (without caldb assistance) processed darks need to be
        manually hardcode swapped based on input science xptime and nd_filter values.
        When MX data is brought into GOA and Caldb, can replace with calls, should
        use caldb argument to figure out which arm is needed instead of here

        Parameters
        ----------
        adinputs: AstroData object(s) for which dark subtraction is to be performed
        dark: (optional) adinput of relevant processed dark
        indiviual: (bool) if True, creates a calib call for each individual science frame.
            If False, groups frames into exposure time and ND_filter and calls one calib per group
            based on first frame in group.

        Returns
        -------
        adinputs with additional image extension of dark subtracted full frame.

        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        adoutputs = []
        if len(adinputs) > 1: # Logic for multiple inputs
            if not individual:  # group frames given by unique sets of xptime and nd_filter
                # saves on redundent caldb requests.
                # Each daytime-made xptime & nd_filter unique processed dark
                # only needs to be called once per science series with that xptime & nd_filter
                exposure_time_list = []
                nd_filter_list = []
                for ad in adinputs:
                    # Create lists of the exposure times and neutral density filters for each image
                    exposure_time_list.append(ad.exposure_time())
                    nd_filter_list.append(ad.filter_orientation()['ND'])

                # Convert lists to numpy arrays for easier indexing
                exposure_time_list = np.array(exposure_time_list)
                nd_filter_list = np.array(nd_filter_list)

                for time in np.unique(exposure_time_list):
                    # Loop over unique exposure times
                    for nd_filter in np.unique(nd_filter_list[exposure_time_list == time]):
                        # Loop over unique ND filters for each exposure time
                        cal_list = []
                        for ad in adinputs:
                            # Create a list of all images with the current exposure time
                            # and ND filter
                            if ad.exposure_time() == time and ad.filter_orientation()['ND'] == nd_filter:
                                cal_list.append(ad)
                        if dark is None:

                            if False:
                                # brainstorming solutions
                                arm_tag = 'BLUE' if 'BLUE' in cal_list[0].tags else 'RED'
                                dark_coeff = _get_calibration_dark_coeff(arm_tag)
                                dark_ad = create_synthetic_dark(dark_coeff, exptime_value=time, nd_value=nd_filter)

                            if 'BLUE' in cal_list[0].tags:
                                # Find the processed dark for the blue arm
                                # TODO: replace hardcoded reference with caldb call to processed dark
                                calib_dark_dir = _get_calibration_dark_path()
                                dark_ad = astrodata.open( str(calib_dark_dir / '20241115T193254Z_DDDDE_b_0300_dark.fits') )
                            elif 'RED' in cal_list[0].tags:
                                # Find the processed dark for the red arm
                                # TODO: replace hardcoded reference with caldb call to processed dark
                                calib_dark_dir = _get_calibration_dark_path()
                                dark_ad = astrodata.open( str(calib_dark_dir / '20241115T193254Z_DDDDE_r_0300_dark.fits') )

                            else:
                                # This condition should never be reached
                                log.warning(f"No dark subtraction will be made to {cal_list[0].filename} "
                                            "group prior to stripe extraction, since no "
                                            "dark was found/specified")

                        for ad_found in cal_list:
                            # Loop over all images in the current exposure time and ND filter group
                            adout = copy.deepcopy(ad_found)
                            log.fullinfo(f"{dark_ad.filename} found as associated dark")
                            # Perform a dark subtraction by using numpy to a do pixel-by-pixel subtraction
                            adout[0].DARK_SUBTRACTED = copy.deepcopy(ad_found)[0].data - copy.deepcopy(dark_ad).data[0]
                            adoutputs.append(adout)
                            if dark_ad:
                                gt.mark_history(ad, primname=self.myself(),
                                                keyword='REDUCTION_DARK',
                                                comment=dark_ad.filename)
        else: # Logic for single input
            for ad in adinputs:
                adout = copy.deepcopy(ad)  # don't group frames and make a dark caldb call for each science frame
                if dark is None:
                    if 'BLUE' in ad.tags:
                        # Find the processed dark for the blue arm
                        # TODO: replace hardcoded reference with caldb call to processed dark
                        calib_dark_dir = _get_calibration_dark_path()
                        dark_ad = astrodata.open( str(calib_dark_dir / '20241115T193254Z_DDDDE_b_0300_dark.fits') )
                    elif 'RED' in ad.tags:
                        # Find the processed dark for the red arm
                        # TODO: replace hardcoded reference with caldb call to processed dark
                        calib_dark_dir = _get_calibration_dark_path()
                        dark_ad = astrodata.open( str(calib_dark_dir / '20241115T193254Z_DDDDE_r_0300_dark.fits') )

                    else:
                        # This condition should never be reached
                        log.warning(f"No dark subtraction will be made to {ad.filename} "
                                    "prior to stripe extraction, since no "
                                    "dark was found/specified")

                log.fullinfo(f"{dark_ad.filename} found as associated dark")
                # Perform a dark subtraction by using numpy to a do pixel-by-pixel subtraction
                adout[0].DARK_SUBTRACTED = copy.deepcopy(ad)[0].data - copy.deepcopy(dark_ad).data[0]
                adoutputs.append(adout)
                if dark_ad:
                    gt.mark_history(ad, primname=self.myself(),
                                    keyword='REDUCTION_DARK',
                                    comment=dark_ad.filename)
            gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adoutputs

    def extractStripes(self, adinputs=None, flat=None,
                       skip_dark=None, slit_height=10,
                       test_extraction=False, individual=False, **params):
        """
        Extracts the stripes from the original 2D spectrum to a sparse array,
        containing only relevant pixels.

        This function marks all relevant pixels for extraction.
        Reinterpreting the flat reference it iterates over all stripes in the
        image and saves a sparse matrix for each stripe.

        TODO: In current format (without caldb assistance) processed flats need to be
        manually hardcoded.
        When MX data is brought into GOA and Caldb, can replace with calls, should
        use caldb argument to figure out which arm is needed instead of here

        Parameters
        ----------
            adinputs
            flat: adinput of relevant processed flat, as processed, will have
                the STRIPES_ID and STRIPES_FIBERS extensions needed
            dark: (optional) adinput of relevant processed dark
            skip_dark: if dark given, which individual fibers dark
                subtraction should be skipped
            slit_height (int): total slit height in px
            test_extraction (bool): used in unit test for this function, saves
            science extraction, flat extraction, and the bpm-extraction in
            FITS-readable format (STRIPES, F_STRIPES, STRIPES_MASK)
            individual: (bool) if False uses one calib call for all frames per arm,
            if True performs a calib call for each frame

        Returns
        -------
            adinputs with sparse matrices added holding the 2D extractions for
            each fiber/order for the science frame, flat frame, and BPM
            (STRIPES, F_STRIPES, STRIPES_MASK)
            if test_extraction==True, the extractions are FITS-readable and not
            sparse matrix format

        """
        if skip_dark is None:
            # skip all dark subtraction by default
            skip_dark = [1, 2, 3, 4, 5]

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        blue_frames = []
        red_frames = []
        for ad in adinputs:
            # create lists of blue and red frames based on tags
            if 'BLUE' in ad.tags:
                blue_frames.append(ad)
            elif 'RED' in ad.tags:
                red_frames.append(ad)

        for ad in blue_frames+red_frames:
            if not ad:
                continue
            try:
                #TODO: I think a lot of this can be replaced with caldb calls - Rohan
                if ad.filename == blue_frames[0].filename:
                    if flat is None:
                        # call one flat for all blue frames, saves redundant caldb requests as new flats are rarely made
                        # time between flats > 1 year
                        # TODO: replace hardcoded reference with caldb call to processed flat
                        calib_flat_dir = _get_calibration_flat_path()
                        flat_ad = astrodata.open( str(calib_flat_dir / '20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits') )
                        log.fullinfo(f"{flat_ad.filename} found as associated flat for all blue frames")
            except IndexError:
                pass
            try:
                if ad.filename == red_frames[0].filename:
                    if flat is None:
                        # call one flat for all red frames, saves redundant caldb requests as new flats are rarely made
                        # time between flats > 1 year
                        # TODO: replace hardcoded reference with caldb call to processed flat
                        calib_flat_dir = _get_calibration_flat_path()
                        flat_ad = astrodata.open( str(calib_flat_dir / '20241114T190714Z_DDDDF_r_0002_DFFFF_flat.fits') )
                        log.fullinfo(f"{flat_ad.filename} found as associated flat for all red frames")
            except IndexError:
                pass #TODO: Implement what will happen here
            if individual:
                # overwrite flat to be used as requested with an individual call for the specific science frame
                if flat is None:
                    if 'BLUE' in ad.tags:
                        # TODO: replace hardcoded reference with caldb call to processed flat
                        calib_flat_dir = _get_calibration_flat_path()
                        flat_ad_i = astrodata.open( str(calib_flat_dir / '20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits') )
                        # '20200911T220106Z_FDDDF_b_0002_FFFFF_flat.fits')
                    elif 'RED' in ad.tags:
                        # TODO: replace hardcoded reference with caldb call to processed flat
                        calib_flat_dir = _get_calibration_flat_path()
                        flat_ad_i = astrodata.open( str(calib_flat_dir / '20241114T190714Z_DDDDF_r_0002_DFFFF_flat.fits') )
                        # '20200911T220106Z_FDDDF_r_0000_FFFFF_flat.fits')
                    else:
                        log.warning(f"No extraction will be made on {ad.filename}, since no "
                                    "flat was found/specified")
                if flat_ad_i:
                    log.fullinfo(f"{flat_ad_i.filename} found as associated flat")
                flat_ad = flat_ad_i
            stripes = {}
            f_stripes = {}
            stripes_masks = {}
            p_id = flat_ad[0].STRIPES_ID

            ad[0].STRIPES_ID = p_id  # record reference to science frame


            # repackage p_id as dict expected format for extraction
            p_id_new = {}
            for i, f in enumerate(flat_ad[0].STRIPES_FIBERS):
                # the table has 6 coefficient per fiber, in ascending order
                p_id_new[f'fiber_{f}'] = dict((colname, p_id[6*i:6*(i+1)][colname].data)
                                        for colname in p_id.colnames)
            p_id = p_id_new


            log.fullinfo("Flat-Identified pixel associations with fiber/order "
                         "found as polynomial info in association "
                         f"with science frame {ad.filename}")

            for f, op in p_id.items():  # extract info into sparse matrices
                adint = copy.deepcopy(ad)
                flatint = copy.deepcopy(flat_ad)
                # dark subtract the frame for the fiber if appropriate
                # each fiber is independent in its need of the dark subtraction
                if int(f[-1]) in skip_dark:
                    log.fullinfo(f'No dark subtracted for fiber {f[-1]}')
                    adint.data[0] = ad.data[0]
                else:
                    #Make a deepcopy of the dark subtracted fiber
                    adint.data[0] = copy.deepcopy(ad)[0].DARK_SUBTRACTED

                log.fullinfo('skipping all fiber dark subtraction is the '
                             'default option')

                for o, p in op.items():
                    # extract the stripe, and the flat stripe and the stripe mask
                    stripe = _extract_single_stripe(
                        adint.data[0], p, slit_height)
                    f_stripe = _extract_single_stripe(
                        flatint.data[0], p, slit_height)
                    s_mask = _extract_single_stripe(
                        np.logical_not(adint.mask[0]).astype(int), p,
                        slit_height)

                    # Update the stripe, flat stripe and stripe mask
                    if f in stripes:
                        stripes[f].update({o: stripe})
                        f_stripes[f].update({o: f_stripe})
                        stripes_masks[f].update({o: s_mask})
                    else:
                        stripes[f] = {o: stripe}
                        f_stripes[f] = {o: f_stripe}
                        stripes_masks[f] = {o: s_mask}
            # store the stripe, flat stripe and stripe mask
            ad[0].STRIPES = stripes
            ad[0].F_STRIPES = f_stripes
            # could find a way to make the flat info retrieve more efficient
            ad[0].STRIPES_MASKS = stripes_masks
            # here mask no longer conforms to BPM, (STRIPES_MASKS == 1 is good)
            # test extraction allows for future reduction direct comparison
            if test_extraction:
                repack_stripes = []
                repack_f_stripes = []
                repack_stripes_masks = []
                test_orders = []
                for ifib in sorted(stripes.keys(), key=lambda x: x.lower()):
                    for iorder in sorted(stripes[ifib].keys(), key=lambda x: x.lower())[:1]:
                        # Convert from sparse matrices to dense matrices
                        repack_stripes.append(stripes[ifib][iorder].todense())
                        repack_f_stripes.append(f_stripes[ifib][iorder].todense())
                        repack_stripes_masks.append(stripes_masks[ifib][iorder].todense())
                        test_orders.append(iorder)

                # Store the stripe, flat stripe and stripe mask as extensions
                ad[0].STRIPES = np.array(repack_stripes)
                ad[0].F_STRIPES = np.array(repack_f_stripes)
                ad[0].STRIPES_MASKS = np.array(repack_stripes_masks)
                ad[0].TEST_ORDERS = np.array(test_orders).astype('int')
            # fix mark history to give full flat and dark name
            gt.mark_history(ad, primname=self.myself(),
                            keyword='REDUCTION_FLAT', comment=flat_ad.filename)
            # if dark_ad:
            #     gt.mark_history(ad, primname=self.myself(),
            #                     keyword='REDUCTION_DARK',
            #                     comment=dark_ad.filename)
            ad.update_filename(suffix=params['suffix'], strip=True)
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def optimalExtraction(self, adinputs=None, opt_extraction=None,
                          back_var=None, full_output=False, penalty=None,
                          s_clip=None, **params):
        """
        Optimal extraction of the 2d echelle spectrum.

        This function performs an optimal extraction of a 2d echelle spectrum.
        A given flat field spectrum is used to generate normalized 'profiles'
        that are used as weighting functions for the spectrum that is going
        to be extracted.
        The algorithm further checks for outliers and rejects them.
        This is to prevent contributions from cosmic hits.

        Parameters
        ----------
            adinputs with STRIPES, F_STRIPES, and STRIPES_MASKS 'extensions' as
            dicts of sparse arrays
            opt_extraction (list): fibers considered for optimal extraction
            back_var (float): manual background variance for frame
            full_output (bool): if True, an additional set of intermediate
                products will be returned / saved
            penalty (float): scaling penalty factor for mismatch correction
                between flat field profile and science spectrum during optimal
                extraction
            s_clip (float): sigma-clipping paramter during optimal extraction

        Returns
        -------
            adinputs with optimal and box extracted orders for each fiber as
            well as uncertainties and the bad pixel mask result from the optimal
            extraction
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # If no fibers are specified, optimal extract fibers 2,3 and 4 (Object fibers in SOOOE frames)
        if opt_extraction is None:
            opt_extraction = [0, 2, 3, 4, 0]
        optimal_reduced_stripes = {}
        box_reduced_stripes = {}
        box_reduced_err = {}
        optimal_reduced_err = {}
        optimal_reduced_2d_arrays = {}
        extracted_bpms = {}

        for ad in adinputs:
            # For each fiber, we need extensions for the reduced orders, the optimal reduced fiber,
            # the error of the optimal reduced fiber, the box reduced fiber, the error of the box
            # reduced fiber, and the bad pixel mask.
            for f in range(1, 6):
                setattr(ad[0], f"REDUCED_ORDERS_FIBER_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"OPTIMAL_REDUCED_FIBER_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"OPTIMAL_REDUCED_ERR_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"BOX_REDUCED_FIBER_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"BOX_REDUCED_ERR_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"BOX_REDUCED_FLAT_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"BPM_FIBER_{f}", np.zeros(shape=[1, 1]))

            # Creating sparse matrices that we end up deleting
            stripes = ad[0].STRIPES
            mask = ad[0].STRIPES_MASKS
            flat_stripes = ad[0].F_STRIPES

            gain = ad.gain()[0][0]  # fix for different channels
            read_noise = ad.read_noise()[0][0]  # Each chip has a different read noise

            for f in stripes.keys():
                if int(f[-1]) in opt_extraction:

                    # if last item of the key is in opt_extraction, we do optimal extraction
                    for o, stripe in stripes[f].items():
                        log.fullinfo(f'Optimum extraction in {f}, order {o}')
                        flux, var, stand_spec, stand_err, fo = _optimal_extraction_single_stripe(
                            stripe, flat_stripes[f][o], gain=gain,
                            read_noise=read_noise, back_var=back_var,
                            mask=mask[f][o], s_clip=s_clip,
                            penalty=penalty, full_output=full_output, log=log)

                        if f in optimal_reduced_stripes:
                            # Update the extensions we created earlier
                            optimal_reduced_stripes[f].update({o: flux})
                            optimal_reduced_err[f].update({o: np.sqrt(var)}) # TODO: When we add the addVAR() method,
                            # we need  to deal with the variance before this line
                            optimal_reduced_2d_arrays[f].update({o: fo})
                            box_reduced_stripes[f].update({o: stand_spec})
                            box_reduced_err[f].update({o: stand_err})
                            extracted_bpms[f].update({o: np.array(np.sum(
                                mask[f][o], axis=0).T).flatten()})
                        else:
                            # We do not have the extensions yet, so create them
                            optimal_reduced_stripes[f] = {o: flux}
                            optimal_reduced_err[f] = {o: var}
                            optimal_reduced_2d_arrays[f] = {o: fo}
                            box_reduced_stripes[f] = {o: stand_spec}
                            box_reduced_err[f] = {o: stand_err}
                            extracted_bpms[f] = {o: np.array(np.sum(
                                mask[f][o], axis=0).T).flatten()}
                else:
                    # if last item of the key is not in opt_extraction, we do box extraction
                    for o, stripe in stripes[f].items():
                        log.fullinfo(f'Only box extraction in {f}, order {o}')
                        stand_spec = _box_extract_single_stripe(
                            stripe, mask[f][o])
                        stand_err = np.sqrt(stand_spec/gain)
                        if f in box_reduced_stripes:
                            # Update the extensions we created earlier
                            box_reduced_stripes[f].update({o: stand_spec})
                            box_reduced_err[f].update({o: stand_err})
                            extracted_bpms[f].update({o: np.array(np.sum(
                                mask[f][o], axis=0).T).flatten()})
                        else:
                            # We do not have the extensions yet, so create them
                            box_reduced_stripes[f] = {o: stand_spec}
                            box_reduced_err[f] = {o: stand_err}
                            extracted_bpms[f] = {o: np.array(np.sum(
                                mask[f][o], axis=0).T).flatten()}
            for f in stripes.keys():
                if f in optimal_reduced_stripes:
                    optimal_reduced_single_fiber = np.array(list(
                        optimal_reduced_stripes[f].values()), dtype=float)
                    optimal_reduced_single_fiber_order_key = np.array(list(
                        optimal_reduced_stripes[f].keys()), dtype=float)
                    optimal_reduced_single_fiber_err = np.array(list(
                        optimal_reduced_err[f].values()), dtype=float)
                    box_reduced_single_fiber = np.array(list(
                        box_reduced_stripes[f].values()), dtype=float)
                    box_reduced_single_err = np.array(list(
                        box_reduced_err[f].values()), dtype=float)
                    bpm_single_fiber = np.array(list(
                        extracted_bpms[f].values()), dtype=int)

                    # Update the extensions based on which fiber we have
                    f_num = int(f[-1]) 
                    setattr(ad[0], f"REDUCED_ORDERS_FIBER_{f_num}", optimal_reduced_single_fiber_order_key)
                    setattr(ad[0], f"OPTIMAL_REDUCED_FIBER_{f_num}", optimal_reduced_single_fiber)
                    setattr(ad[0], f"OPTIMAL_REDUCED_ERR_{f_num}", optimal_reduced_single_fiber_err)
                    setattr(ad[0], f"BOX_REDUCED_FIBER_{f_num}", box_reduced_single_fiber)
                    setattr(ad[0], f"BOX_REDUCED_ERR_{f_num}", box_reduced_single_err)
                    # setattr(ad[0], f"BOX_REDUCED_FLAT_{f_num}", np.zeros(shape=[1, 1]))
                    setattr(ad[0], f"BPM_FIBER_{f_num}", bpm_single_fiber)
                else:
                    # Dealing with the case that we have no optimal extraction, so the reduced order
                    # is the box extraction as opposed to the optimal extraction.
                    box_reduced_single_fiber_order_key = np.array(list(
                        box_reduced_stripes[f].keys()), dtype=float)
                    box_reduced_single_fiber = np.array(list(
                        box_reduced_stripes[f].values()), dtype=float)
                    box_reduced_single_err = np.array(list(
                        box_reduced_err[f].values()), dtype=float)
                    bpm_single_fiber = np.array(list(
                        extracted_bpms[f].values()), dtype=int)

                    # Update the extensions based on which fiber we have.
                    f_num = int(f[-1]) 
                    setattr(ad[0], f"REDUCED_ORDERS_FIBER_{f_num}", box_reduced_single_fiber_order_key)
                    setattr(ad[0], f"BOX_REDUCED_FIBER_{f_num}", box_reduced_single_fiber)
                    setattr(ad[0], f"BOX_REDUCED_ERR_{f_num}", box_reduced_single_err)
                    # setattr(ad[0], f"BOX_REDUCED_FLAT_{f_num}", box_reduced_flat)
                    setattr(ad[0], f"BPM_FIBER_{f_num}", bpm_single_fiber)

            # Delete the sparse matrices from the ad object
            del ad[0].STRIPES
            del ad[0].F_STRIPES
            del ad[0].STRIPES_MASKS

            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=False)
            log.fullinfo(f"frame {ad.filename} extracted")
        return adinputs

    def boxExtraction(self, adinputs, **params):
        """
        This primitive performs box extraction on a 2d echelle spectrum.
        Utilized in the dynamic and static wavelength calibration recipe as it
        is quicker than relying on optimal extraction.

        Parameters
        ----------
        adinputs with STRIPES, F_STRIPES, and STRIPES_MASKS 'extensions' as
            dicts of sparse arrays

        Returns
        -------
        adinputs with box extracted orders for each fiber as
        well as uncertainties calculated during the box extraction
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        box_reduced_stripes = {}
        box_reduced_err = {}
        box_reduced_flats = {}
        extracted_bpms = {}

        for ad in adinputs:
            # For each fiber, we need extensions for the reduced orders,
            # the box reduced fiber, the error of the boxreduced fiber,
            # and the bad pixel mask.
            for f in range(1, 6):
                setattr(ad[0], f"REDUCED_ORDERS_FIBER_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"BOX_REDUCED_FIBER_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"BOX_REDUCED_ERR_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"BOX_REDUCED_FLAT_{f}", np.zeros(shape=[1, 1]))
                setattr(ad[0], f"BPM_FIBER_{f}", np.zeros(shape=[1, 1]))

            # Creating sparse matrices
            stripes = ad[0].STRIPES
            mask = ad[0].STRIPES_MASKS
            flat_stripes = ad[0].F_STRIPES

            gain = ad.gain()[0][0]  # fix for different channels
            read_noise = ad.read_noise()[0][0]  # Each chip has a different read noise

            for f in stripes:
                for o, stripe in stripes[f].items():
                    log.fullinfo(f'Box extraction in {f}, order {o}')
                    
                    stand_spec = _box_extract_single_stripe(stripe, mask[f][o])
                    stand_err = np.sqrt(stand_spec/gain)
                    stand_flat = _box_extract_single_stripe(flat_stripes[f][o], mask[f][o])
                    if f in box_reduced_stripes:
                        # Update the extensions we created earlier
                        box_reduced_stripes[f].update({o: stand_spec})
                        box_reduced_err[f].update({o: stand_err})
                        box_reduced_flats[f].update({o: stand_flat})
                        extracted_bpms[f] = {o: np.array(np.sum(mask[f][o], axis=0).T).flatten()}
                    else:
                        # We do not have the extensions yet, so create them
                        box_reduced_stripes[f] = {o: stand_spec}
                        box_reduced_err[f] = {o: stand_err}
                        box_reduced_flats[f] = {o: stand_flat}
                        extracted_bpms[f] = {o: np.array(np.sum(mask[f][o], axis=0).T).flatten()}

                box_reduced_single_fiber_order_key = np.array(list(
                    box_reduced_stripes[f].keys()), dtype=float)
                box_reduced_single_fiber = np.array(list(
                        box_reduced_stripes[f].values()), dtype=float)
                box_reduced_single_err = np.array(list(
                        box_reduced_err[f].values()), dtype=float)
                box_reduced_flat = np.array(list(
                        box_reduced_flats[f].values()), dtype=float)
                bpm_single_fiber = np.array(list(extracted_bpms[f].values()), dtype=int)

                # Update the extensions based on which fiber we have
                f_num = int(f[-1])  # get the fiber number from the key    
                setattr(ad[0], f"REDUCED_ORDERS_FIBER_{f_num}", box_reduced_single_fiber_order_key)
                setattr(ad[0], f"BOX_REDUCED_FIBER_{f_num}", box_reduced_single_fiber)
                setattr(ad[0], f"BOX_REDUCED_ERR_{f_num}", box_reduced_single_err)
                setattr(ad[0], f"BOX_REDUCED_FLAT_{f_num}", box_reduced_flat)
                setattr(ad[0], f"BPM_FIBER_{f_num}", bpm_single_fiber)

            del ad[0].STRIPES, ad[0].F_STRIPES, ad[0].STRIPES_MASKS
            ad.update_filename(suffix=params["suffix"], strip=False)
            log.fullinfo(f"frame {ad.filename} extracted")
        return adinputs


##############################################################################
# Below are the helper functions for the primitives in this module           #
##############################################################################

def create_synthetic_dark(ad_coeff, exptime_value=None, nd_value=None):
    """
    Create a synthetic dark frame from coefficients file.
    
    This function generates a synthetic dark frame using pre-computed coefficients
    and the log-linear relationship: dark = z1 + z0 * log10(exptime * factor)
    
    Parameters
    ----------
    ad_coeff : AstroData
        AstroData object containing coefficient data with extensions:
        - COEFF_Z0: slope coefficients for each pixel
        - COEFF_Z1: intercept coefficients for each pixel  
        - LOGEXPTIME: table with exposure times and ND filter data
    exptime_value : float, optional
        Target exposure time in seconds. If None, must provide nd_value.
    nd_value : float, optional
        Target ND filter position. If None, will be calculated from exptime_value.
        
    Returns
    -------
    synthetic_dark : numpy.ndarray
        2D array containing the synthetic dark frame
    actual_exptime : float
        The actual exposure time used (may differ from input due to ND corrections)
    actual_nd : float
        The actual ND filter position used
    factor : float
        The correction factor applied (1.0 if no correction needed)
        
    Raises
    ------
    ValueError
        If neither exptime_value nor nd_value is provided, or if required 
        extensions are missing from ad_coeff
    """
    
    # Validate inputs
    if exptime_value is None and nd_value is None:
        raise ValueError("Either exptime_value or nd_value must be provided")
    
    # Check for required extensions
    required_extensions = ['COEFF_Z0', 'COEFF_Z1', 'LOGEXPTIME']
    for ext_name in required_extensions:
        if not hasattr(ad_coeff[0], ext_name):
            raise ValueError(f"Required extension {ext_name} not found")
    
    # Extract coefficient arrays
    z0 = ad_coeff[0].COEFF_Z0
    z1 = ad_coeff[0].COEFF_Z1
    logexptime_table = ad_coeff[0].LOGEXPTIME
    
    # Extract calibration data
    logexptimes = np.array(logexptime_table['logexptime'])
    exptimes = np.array(logexptime_table['exptime'])
    
    # Check if ND filter data is available
    has_nd_data = 'ndfilter' in logexptime_table.colnames
    if has_nd_data:
        ndfilters = np.array(logexptime_table['ndfilter'])
    else:
        ndfilters = np.zeros_like(exptimes)  # Default to zero if no ND data
    
    # Initialize variables
    factor = 1.0
    actual_exptime = exptime_value
    actual_nd = nd_value
    
    # Case 1: Only exposure time provided
    if exptime_value is not None and nd_value is None:
        actual_exptime = exptime_value
        
        if has_nd_data and len(np.unique(ndfilters)) > 1:
            # Calculate expected ND position from exposure time
            z_nd_logt = np.polyfit(logexptimes, ndfilters, 1)
            f_nd_logt = np.poly1d(z_nd_logt)
            actual_nd = f_nd_logt(np.log10(exptime_value))
        else:
            actual_nd = ndfilters[0] if len(ndfilters) > 0 else 0.0
    
    # Case 2: Only ND value provided  
    elif nd_value is not None and exptime_value is None:
        actual_nd = nd_value
        
        if has_nd_data and len(np.unique(ndfilters)) > 1:
            # Calculate exposure time from ND position
            z_logt_nd = np.polyfit(ndfilters, logexptimes, 1)
            f_logt_nd = np.poly1d(z_logt_nd)
            actual_exptime = 10**(f_logt_nd(nd_value))
        else:
            # Use median exposure time as default
            actual_exptime = np.median(exptimes)
    
    # Case 3: Both exposure time and ND value provided
    elif exptime_value is not None and nd_value is not None:
        actual_exptime = exptime_value
        actual_nd = nd_value
        
        # Check for ND filter mismatch and apply correction if needed
        if has_nd_data and len(np.unique(ndfilters)) > 1:
            z_nd_logt = np.polyfit(logexptimes, ndfilters, 1)
            f_nd_logt = np.poly1d(z_nd_logt)
            z_logt_nd = np.polyfit(ndfilters, logexptimes, 1) 
            f_logt_nd = np.poly1d(z_logt_nd)
            
            expected_nd = f_nd_logt(np.log10(exptime_value))
            nd_difference = abs(actual_nd - expected_nd)
            
            # Apply correction factor if ND mismatch is significant (>0.2)
            if nd_difference > 0.2:
                logt_nominal = f_logt_nd(actual_nd)
                factor = 10**(np.log10(exptime_value) - logt_nominal)
    
    # Calculate synthetic dark frame
    # Formula: dark = z1 + z0 * log10(exptime * factor)
    effective_logexptime = np.log10(actual_exptime * factor)
    synthetic_dark = z1 + z0 * effective_logexptime
    
    return synthetic_dark.astype(np.float32), actual_exptime, actual_nd, factor


def _extract_single_stripe(data=None, polynomials=None, slit_height=10):
    """
    Extracts single stripe from 2d image.

    This function returns a sparse matrix containing all relevant pixels for
    a single stripe for a given polynomial p and a given slit height.

    Parameters
    ----------
        polynomials (np.ndarray): polynomial coefficients
        img (np.ndarray): 2d echelle spectrum
        slit_height (int): total slit height in pixel to be extracted

    Returns
    -------
        scipy.sparse.csc_matrix: extracted spectrum

    """
    ny, nx = data.shape
    if isinstance(slit_height, np.ndarray):
        slit_height = int(slit_height[0])

    xx = np.arange(nx)
    y = np.poly1d(polynomials)(xx)

    # Create matrix containing all indices for a given slit height.
    # For the y matrix, the matrix is values from - slit height to + slit height,
    # repeating nx times.  For the x matrix, it is values 0 to nx repeated 2* slit height.
    # Both have dimensions (2*slit height) x nx.

    slit_indices_y = np.arange(-slit_height, slit_height
                                ).repeat(nx).reshape((2 * slit_height, nx))
    slit_indices_x = np.tile(np.arange(nx), 2 * slit_height
                                ).reshape((2 * slit_height, nx))

    indices = np.rint(slit_indices_y + y).astype(int)
    valid_indices = np.logical_and(indices < ny, indices > 0)

    # Create sparse matrix of dimensions ny x nx, with row indices given by indices[valid_indices]
    # and column indices sit_indices_x[valid_indices].

    mat = sparse.coo_matrix((data[indices[valid_indices],
        slit_indices_x[valid_indices]],(indices[valid_indices],
        slit_indices_x[valid_indices])), shape=(ny, nx))
    return mat.tocsc()

@numba.njit(cache=True)
def reject_step(back_var, data_var, debug_level, diff_aver, diff_save, flux, gain, good_disp, mask, profile, read_noise,
            reject_tracker, s_clip, stand_spec, stripe, var):
    for h in good_disp:
        expected = profile[:, h] * mask[:, h] * stand_spec[h]  # flat column with mask scaled to data total
        actual = stripe[:, h] * mask[:, h]  # actual column data with mask
        diff = actual - expected

        data_var[:, h] = np.abs(stand_spec[h] * profile[:, h]) / gain + back_var[:, h] + read_noise
        noise_rev = 1 / np.sqrt(data_var[:, h])
        diff_save[:, h] = diff * noise_rev
        reject_index = np.nonzero(((np.abs(diff) - np.abs(diff_aver[:, h])) * noise_rev) >= s_clip)[0]

        while len(reject_index) > 0:
            worst = np.argmax((np.abs(diff) - np.abs(diff_aver[:, h])) * noise_rev)
            reject_tracker[h] = reject_tracker[h] + 1
            # if debug_level >= 3:
            #     logger.info(f'Outlier found in column {h}, pixel {worst}')
            #     fig, ax = plt.subplots(2, 1)
            #     actual_plot = actual.copy()
            #     expected_plot = expected.copy()
            #     actual_plot[actual_plot == 0] = np.nan
            #     expected_plot[expected_plot == 0] = np.nan
            #     ax[0].plot(actual_plot, 'r')
            #     ax[0].plot(expected_plot, 'b')
            #     ax[1].plot(np.abs(diff), 'r')
            #     ax[1].plot(np.abs(diff) - np.abs(diff_aver[:, h]), 'g')
            #     ax[1].plot(np.sqrt(data_var[:, h]) * s_clip, 'b')
            #     plt.show()

            mask[worst, h] = 0

            denom = np.sum(profile[:, h] * profile[:, h] * mask[:, h] / data_var[:, h])

            stand_spec[h] = np.sum(profile[:, h] * mask[:, h] * stripe[:, h] / data_var[:, h]) / denom

            expected = profile[:, h] * mask[:, h] * stand_spec[h]
            actual = stripe[:, h] * mask[:, h]
            diff = actual - expected

            data_var[:, h] = np.abs(stand_spec[h] * profile[:, h]) / gain + back_var[:, h] + read_noise
            noise_rev = 1 / np.sqrt(data_var[:, h])
            # diff_save[:,h] = diff*noise_rev
            reject_index = np.nonzero(((np.abs(diff) - np.abs(diff_aver[:, h])) * noise_rev) >= s_clip)[0]

            if np.count_nonzero(actual[3:-5]) < len(profile[3:-5, h]) / 2.:
                # reject_index = np.array([], dtype=int)
                reject_index = np.empty((0), dtype=np.int_)
                mask[:, h] = 0
                flux[h] = 0
            # logger.warning(f'Too many bad pixels in column {h}, reject column')
        if np.count_nonzero(mask[:, h]) > 0:
            denom = np.sum(profile[:, h] * profile[:, h] * mask[:, h] / data_var[:, h])
            flux[h] = np.sum(profile[:, h] * mask[:, h] * stripe[:, h] / data_var[:, h]) / denom
            var[h] = np.sum(profile[:, h] * mask[:, h]) / denom

#@staticmethod
def rectify_and_reshape(back_var, box_extracted_flat_stripe, data, flat_stripe, mask, new_back_var, new_mask, profile,
                    stripe):
    """
    This function takes the stripe and mask and reshapes them to be used in the optimal extraction algorithm.
    It also creates the profile and background variance arrays.

    Parameters
    ----------
    """
    first_row_bigger_than_0 = (stripe!=0).argmax(axis=0)  # returns 0 for rows that contain only 0s
    slit_height = data.shape[0]
    data[:,:] = np.vstack([stripe[fb:fb + slit_height, i] for i, fb in enumerate(first_row_bigger_than_0)]).T
    new_mask[:,:] = np.vstack([mask[fb:fb + slit_height, i] for i, fb in enumerate(first_row_bigger_than_0)]).T
    new_back_var[:,:] = np.vstack([back_var[fb:fb + slit_height, i] for i, fb in enumerate(first_row_bigger_than_0)]).T
    profile[:,:] = np.vstack([flat_stripe[fb:fb + slit_height, i]/box_extracted_flat_stripe[i] for i, fb in enumerate(first_row_bigger_than_0)]).T

    print(slit_height)

def _optimal_extraction_single_stripe(stripe, flat_stripe, gain=1, read_noise=1.23, back_var=None, mask=None,
                                    debug_level=0, full_output=False, s_clip=5.0, penalty=1.0, log=None):
    """
    Performs optimal extraction of a single stripe.
    Based on the algorithm described by Horne et al. 1986, PASP, 98, 609.

    Args:
        stripe (scipy.sparse.spmatrix): science frame stripe to be extracted
        flat_stripe (scipy.sparse.spmatrix): flat frame stripe to be used as profile
        gain (float): detector gain factor (conversion photons -> DN,  given in e-/DN )
        read_noise (float): typical detector read noise given as the variance in DN
        back_var (np.ndarray): background variance (from scattered light model oder dark, otherwhise 0)
        mask (np.ndarray): bad pixel mask. Note: the mask will be modified by iterative algorithm that looks for
        outliers.
        debug_level (int): debug level
        full_output (bool): if True, returns all intermediate results for debugging/testing
        s_clip (float): sigma clipping value for optimal extraction
        penalty (float): scaling factor for global per-order profile mismatch correction. Set 0 for no correction

    Returns
    -------
        tuple(np.ndarray, np.ndarray, dict): (optimal extracted spectrum, box extracted spectrum, dict of additional
        intermediate results if full_output was True)

    """
    #log = self.log
    if mask is None:
        mask       = stripe.copy()
        mask[:, :] = 1
    else:
        # copy mask, because it will be modified
        mask     = mask.copy()
        #back_var =  back_var.copy()
        back_var = stripe.copy()
        back_var[:, :] = 0  # back_var.copy()

    # box extracted spectrum
    stand_spec0 = _box_extract_single_stripe(stripe, mask)  # direct box extraction,
    stand_spec = stand_spec0.copy()
    stand_err = np.sqrt(stand_spec / gain)  # error of box extracted spectrum

    # flat data
    box_extracted_flat_stripe = _box_extract_single_stripe(flat_stripe, mask)  # spatial sums for flat along disp

    print(f"stripe shape: {stripe.shape}, flat stripe shape: {flat_stripe.shape}, mask shape: {mask.shape}"
          f", back_var shape: {back_var.shape}, box_extracted_flat_stripe shape: {box_extracted_flat_stripe.shape}")
    # print the types
    print(f"stripe type: {type(stripe)}, flat stripe type: {type(flat_stripe)}, mask type: {type(mask)}"
            f", back_var type: {type(back_var)}, box_extracted_flat_stripe type: {type(box_extracted_flat_stripe)}")
    # cut stripe sparse matrix into numpy array
    # find the spatial columns utilized along entire stripe (greater than slit height because of stripe path)
    sparse_vcols = np.array(~np.all(stripe.todense() == 0, axis=1)).reshape(-1)
    stripe       = np.array(stripe.todense()[sparse_vcols])  # strip stripe to the inclusive nonzero rows

    # mask         = mask[sparse_vcols]  # strip mask similarly
    # back_var     = back_var[sparse_vcols] #strip background variance map similarly
        
    mask         = np.array(mask.todense()[sparse_vcols])  # strip mask similarly
    back_var     = np.array(back_var.todense()[sparse_vcols]) #strip background variance map similarly
    
    flat_stripe  = np.array(flat_stripe.todense()[sparse_vcols])
    sparse_vrows = np.count_nonzero((stripe != 0).T[1500])  # use ~middle column slit height as slit height pass
    data         = np.zeros((sparse_vrows, stripe.shape[1]))
    diff_save    = data.copy()
    new_mask     = data.copy()  # create actual limit numpy arrays
    new_back_var = data.copy()  # create actual limit numpy arrays
    profile      = data.copy()
    
    print(f"stripe shape: {stripe.shape}, flat stripe shape: {flat_stripe.shape}, mask shape: {mask.shape}"
          f", back_var shape: {back_var.shape}")
    # print the types
    print(f"stripe type: {type(stripe)}, flat stripe type: {type(flat_stripe)}, mask type: {type(mask)}"
            f", back_var type: {type(back_var)}")
    #import ipdb; ipdb.set_trace()

    # legacy code that could be optimized
    for i in np.arange(stripe.shape[1]):
        if np.nonzero((stripe != 0).T[i])[0].shape[0] == data.shape[0]:  # if column is slit height (not edge of chip)
            if box_extracted_flat_stripe[i] > 1E-12:
                data[:, i]        = stripe[np.nonzero((stripe != 0).T[i])[0], i]  # write data
                new_mask[:, i]    = mask[np.nonzero((stripe != 0).T[i])[0], i]
                new_back_var[:,i] =  back_var[np.nonzero((stripe != 0).T[i])[0], i]
                profile[:, i]     = flat_stripe[np.nonzero((stripe != 0).T[i])[0], i] / box_extracted_flat_stripe[i]
        else:
            new_mask[:, i]     = 0  # could be optimized
            new_back_var[:, i] = 0  # could be optimized
    
    # TODO: possible optimization, refacter after regression tests pass
    # rectify_and_reshape(back_var, box_extracted_flat_stripe, data, flat_stripe, mask, new_back_var, new_mask, profile, stripe)

    new_mask[:, stand_spec0 < 1E-12] = 0  # if sum is less than zero, whole column is cancelled for this stripe
    mask = new_mask.copy()
    back_var = new_back_var.copy() * mask
    stripe = data  # return variables 'mask' and 'stripe' to the naming conventions
    data_var = abs(stripe.copy())/gain + back_var + read_noise

    # final output
    flux = np.zeros(len(stand_spec0))
    var = flux.copy()

    # Calculate a first guess of the difference between stripe and scaled flat_stripe
    expected = profile * mask * stand_spec
    actual = stripe * mask
    if len(actual) !=len(expected):
        log.warning("Hit flat/science mismatch")
    diff = actual - expected

    # Calculate the median of the difference along the order. This represents the 'global' mismatch between flat and science profile
    # and helps correct for the 'drift' problem in x-dispersion.
    diff_aver = penalty * np.abs(median_filter(diff,size=(1,201)))

    # # avoid already caught bad pixels (whole column is zero in bpm or stripe is too small)
    good_disp = np.nonzero(np.array(~np.any(mask == 0, axis=0)))[0]
    reject_tracker = np.zeros_like(stand_spec0)

    # legacy code that could be optimized
    for h in good_disp:
        expected = profile[:, h] * mask[:, h] * stand_spec[h]  # flat column with mask scaled to data total
        actual = stripe[:, h] * mask[:, h]  # actual column data with mask
        diff = actual - expected

        data_var[:, h] = abs(stand_spec[h] * profile[:, h]) / gain + back_var[:, h] + read_noise
        noise_rev = 1 / np.sqrt(data_var[:, h])
        diff_save[:, h] = diff * noise_rev
        reject_index = np.nonzero(((np.abs(diff) - np.abs(diff_aver[:, h])) * noise_rev) >= s_clip)[0]

        while len(reject_index) > 0:

            worst = np.argmax((np.abs(diff) - np.abs(diff_aver[:, h])) * noise_rev)
            reject_tracker[h] = reject_tracker[h] + 1

            mask[worst, h] = 0

            denom = np.sum(profile[:, h] * profile[:, h] * mask[:, h] / data_var[:, h])

            stand_spec[h] = np.sum(profile[:, h] * mask[:, h] * stripe[:, h] / data_var[:, h]) / denom

            expected = profile[:, h] * mask[:, h] * stand_spec[h]
            actual = stripe[:, h] * mask[:, h]
            diff = actual - expected

            data_var[:, h] = abs(stand_spec[h] * profile[:, h]) / gain + back_var[:, h] + read_noise
            noise_rev = 1 / np.sqrt(data_var[:, h])
            # diff_save[:,h] = diff*noise_rev
            reject_index = np.nonzero(((np.abs(diff) - np.abs(diff_aver[:, h])) * noise_rev) >= s_clip)[0]

            if np.count_nonzero(actual[3:-5]) < len(profile[3:-5, h]) /2.:
                reject_index = np.array([])
                mask[:, h] = 0
                flux[h] = 0
                log.warning(f'Too many bad pixels in column {h}, reject column')
        if np.count_nonzero(mask[:, h]) > 0:
            denom = np.sum(profile[:, h] * profile[:, h] * mask[:, h] / data_var[:, h])
            flux[h] = np.sum(profile[:, h] * mask[:, h] * stripe[:, h] / data_var[:, h]) / denom
            var[h] = np.sum(profile[:, h] * mask[:, h]) / denom

    # TODO: possible optimization, refacter after regression tests pass
    # reject_step(back_var, data_var, debug_level, diff_aver, diff_save, flux, gain, good_disp, mask, profile, read_noise,
    #         reject_tracker, s_clip, stand_spec, stripe, var)

    flux[flux == 0] = np.nan
    var[var == 0] = np.nan

    stand_spec0[stand_spec0 == 0] = np.nan
    stand_err[stand_err == 0] = np.nan

    total_count = np.sum(reject_tracker > 0)
    substantial_count = np.sum(np.abs(stand_spec0 - stand_spec)[reject_tracker > 0] > 0.005 * (stand_spec0[reject_tracker > 0]))
    irrelevant_count = total_count - substantial_count

    if total_count > 0:
        log.fullinfo(f'Rejected {np.sum(reject_tracker):.0f} pixels in {total_count} columns during optimal extraction.')
        irrelevant_count_percentage = irrelevant_count/total_count*100
        if irrelevant_count_percentage > 30:
            log.warning(f'Rejections with flux changes < 0.5%: {irrelevant_count} ({irrelevant_count_percentage:.0f}%)')
        else:
            log.fullinfo(f'Rejections with flux changes < 0.5%: {irrelevant_count} ({irrelevant_count_percentage:.0f}%)')
    else:
        log.fullinfo('No rejections')

    # return all intermediate results for debugging/testing
    if full_output:
        return flux, var, stand_spec0, stand_err, {'noise': var, 'acceptancemask': mask, "stripe": stripe,
                                    "initial_sigma": diff_save}  # {'noise': var, 'profile': profile, 'rejectionmask': sparse.csr_matrix(mask),
    # 'expected': expected, 'actual': actual}
    return flux, var, stand_spec0, stand_err, {'noise': var}

#@staticmethod
def _box_extract_single_stripe(stripe=None, mask=None):
    """
    Box extraction of a single stripe.

    Parameters
    ----------
        stripe (sparse.matrix): stripe to be extracted
        mask (np.ndarray): binary pixel mask. Where mask==1: values will
        contribute. Same shape as stripe.toarray()

    Returns
    -------
        np.ndarray: box extracted flux
    """
    if mask is None:
        mask = stripe.copy()
        mask.data[:] = 1.
    return np.array(np.sum(stripe.multiply(mask), axis=0).T).flatten()
