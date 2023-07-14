#
#                                                                       DRAGONS
#
#                                                 primitives_maroonx_echelle.py
# ------------------------------------------------------------------------------

# from geminidr.gemini.lookups import DQ_definitions as DQ
import copy
import astrodata
import numpy as np
from scipy.ndimage import median_filter
import scipy.sparse as sparse

from .primitives_maroonx import MAROONX
from . import parameters_maroonx_echelle
from recipe_system.utils.decorators import parameter_override
from gempy.gemini import gemini_tools as gt
from geminidr.core import Spect
# ------------------------------------------------------------------------------

@parameter_override
class MAROONXEchelle(MAROONX, Spect):
    """
    This class contains primitives that applies to all MAROON-X echelle
    data.
    """

    tagset = {'GEMINI', 'MAROONX', 'ECHELLE', 'SPECT'}

    def __init__(self, adinputs, **kwargs):
        super(MAROONXEchelle, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'maroonxdr.maroonx.lookups'
        self._param_update(parameters_maroonx_echelle)

    def darkSubtraction(self, adinputs=None, dark=None, individual=False,
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
                # saves on redundent caldb requests. Each daytime-made xptime & nd_filter unique processed dark
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
                            # Create a list of all images with the current exposure time and ND filter
                            if ad.exposure_time() == time and ad.filter_orientation()['ND'] == nd_filter:
                                cal_list.append(ad)
                        if dark is None:
                            if 'BLUE' in cal_list[0].tags:
                                # Find the processed dark for the blue arm
                                # TODO: replace hardcoded reference with caldb call to processed dark
                                dark_ad = astrodata.open('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/calibrations/processed_dark/' +
                                                         '20220808T163524Z_DDDDE_b_0300_dark.fits')
                            elif 'RED' in cal_list[0].tags:
                                # Find the processed dark for the red arm
                                # TODO: replace hardcoded reference with caldb call to processed dark
                                dark_ad = astrodata.open('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/calibrations/processed_dark/' +
                                                         '20220808T163524Z_DDDDE_r_0300_dark.fits')

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
                        dark_ad = astrodata.open('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/calibrations/processed_dark/' +
                                                 '20220808T163524Z_DDDDE_b_0300_dark.fits')
                    elif 'RED' in ad.tags:
                        # Find the processed dark for the red arm
                        # TODO: replace hardcoded reference with caldb call to processed dark
                        dark_ad = astrodata.open('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/calibrations/processed_dark/' +
                                                 '20220808T163524Z_DDDDE_r_0300_dark.fits')

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
                        flat_ad = astrodata.open('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/calibrations/processed_flat/' +
                                                 '20220725T164012Z_FDDDF_b_0007_FFFFF_flat.fits')
                        log.fullinfo(f"{flat_ad.filename} found as associated flat for all blue frames")
            except IndexError:
                pass
            try:
                if ad.filename == red_frames[0].filename:
                    if flat is None:
                        # call one flat for all red frames, saves redundant caldb requests as new flats are rarely made
                        # time between flats > 1 year
                        # TODO: replace hardcoded reference with caldb call to processed flat
                        flat_ad = astrodata.open('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/calibrations/processed_flat/' +
                                                   '20220725T164012Z_FDDDF_r_0001_FFFFF_flat.fits')
                        log.fullinfo(f"{flat_ad.filename} found as associated flat for all red frames")
            except IndexError:
                pass
            if individual:
                # overwrite flat to be used as requested with an individual call for the specific science frame
                if flat is None:
                    if 'BLUE' in ad.tags:
                        # TODO: replace hardcoded reference with caldb call to processed flat
                        flat_ad_i = astrodata.open('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/calibrations/processed_flat/'+
                        '20220725T164012Z_FDDDF_b_0007_FFFFF_flat.fits')
                        # '20200911T220106Z_FDDDF_b_0002_FFFFF_flat.fits')
                    elif 'RED' in ad.tags:
                        # TODO: replace hardcoded reference with caldb call to processed flat
                        flat_ad_i = astrodata.open('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/calibrations/processed_flat/'+
                        '20220725T164012Z_FDDDF_r_0001_FFFFF_flat.fits')
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
            fiber_check = flat_ad[0].STRIPES_FIBERS
            if len(p_id) != 30 or len(fiber_check) != 5:
                log.error(f"No extraction will be made on {ad.filename}, since "
                            "flat is not fully processed")

            ad[0].STRIPES_ID = p_id  # record reference to science frame
            #repackage p_id as dict expected format for extraction
            p_id = {'fiber_1': dict((colname, p_id[0:6][colname].data)
                                    for colname in p_id.colnames),
                        'fiber_2': dict((colname, p_id[6:12][colname].data)
                                        for colname in p_id.colnames),
                        'fiber_3': dict((colname, p_id[12:18][colname].data)
                                        for colname in p_id.colnames),
                        'fiber_4': dict((colname, p_id[18:24][colname].data)
                                        for colname in p_id.colnames),
                        'fiber_5': dict((colname, p_id[24:][colname].data)
                                        for colname in p_id.colnames)}
            
            log.fullinfo("Flat-Identified pixel associations with fiber/order "
                         "found as polynomial info in association "
                         f"with science frame {ad.filename}")
            
            for f in p_id.keys():  # extract info into sparse matrices
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
                
                for o, p in p_id[f].items():
                    # extract the stripe, and the flat stripe and the stripe mask
                    stripe = self._extract_single_stripe(
                        adint.data[0], p, slit_height)
                    f_stripe = self._extract_single_stripe(
                        flatint.data[0], p, slit_height)
                    s_mask = self._extract_single_stripe(
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
        if opt_extraction is None:
            opt_extraction = [0, 2, 3, 4, 0]
        optimal_reduced_stripes = {}
        box_reduced_stripes = {}
        box_reduced_err = {}
        optimal_reduced_err = {}
        optimal_reduced_2d_arrays = {}
        extracted_bpms = {}
        for ad in adinputs:
            """
            For each fiber, we need extensions for the reduced orders, the optimal reduced fiber,
            the error of the optimal reduced fiber, the box reduced fiber, the error of the box
            reduced fiber, and the bad pixel mask.
            """
            # Fiber 1
            ad[0].REDUCED_ORDERS_FIBER_1 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_FIBER_1 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_ERR_1 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_FIBER_1 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_ERR_1 = np.zeros(shape=[1, 1])
            ad[0].BPM_FIBER_1 = np.zeros(shape=[1, 1])
            # Fiber 2
            ad[0].REDUCED_ORDERS_FIBER_2 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_FIBER_2 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_ERR_2 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_FIBER_2 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_ERR_2 = np.zeros(shape=[1, 1])
            ad[0].BPM_FIBER_2 = np.zeros(shape=[1, 1])
            # Fiber 3
            ad[0].REDUCED_ORDERS_FIBER_3 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_FIBER_3 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_ERR_3 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_FIBER_3 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_ERR_3 = np.zeros(shape=[1, 1])
            ad[0].BPM_FIBER_3 = np.zeros(shape=[1, 1])
            # Fiber 4
            ad[0].REDUCED_ORDERS_FIBER_4 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_FIBER_4 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_ERR_4 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_FIBER_4 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_ERR_4 = np.zeros(shape=[1, 1])
            ad[0].BPM_FIBER_4 = np.zeros(shape=[1, 1])
            # Fiber 5
            ad[0].REDUCED_ORDERS_FIBER_5 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_FIBER_5 = np.zeros(shape=[1, 1])
            ad[0].OPTIMAL_REDUCED_ERR_5 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_FIBER_5 = np.zeros(shape=[1, 1])
            ad[0].BOX_REDUCED_ERR_5 = np.zeros(shape=[1, 1])
            ad[0].BPM_FIBER_5 = np.zeros(shape=[1, 1])
            
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
                        flux, var, stand_spec, stand_err, fo = self._optimal_extraction_single_stripe(
                        stripe, flat_stripes[f][o], gain=gain,
                        read_noise=read_noise, back_var=back_var,
                        mask=mask[f][o], s_clip=s_clip,
                        penalty=penalty, full_output=full_output)

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
                        stand_spec = self._box_extract_single_stripe(
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
                if f in optimal_reduced_stripes.keys():
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
                    
                    # Update the extensions depending on the fiber.  A lot of code repetition here, there may be a
                    # cleaner way to do this.
                    if f == 'fiber_1':
                        ad[0].REDUCED_ORDERS_FIBER_1 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_1 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_1 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_1 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_1 = box_reduced_single_err
                        ad[0].BPM_FIBER_1 = bpm_single_fiber
                    if f == 'fiber_2':
                        ad[0].REDUCED_ORDERS_FIBER_2 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_2 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_2 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_2 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_2 = box_reduced_single_err
                        ad[0].BPM_FIBER_2 = bpm_single_fiber
                    if f == 'fiber_3':
                        ad[0].REDUCED_ORDERS_FIBER_3 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_3 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_3 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_3 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_3 = box_reduced_single_err
                        ad[0].BPM_FIBER_3 = bpm_single_fiber
                    if f == 'fiber_4':
                        ad[0].REDUCED_ORDERS_FIBER_4 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_4 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_4 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_4 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_4 = box_reduced_single_err
                        ad[0].BPM_FIBER_4 = bpm_single_fiber
                    if f == 'fiber_5':
                        ad[0].REDUCED_ORDERS_FIBER_5 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_5 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_5 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_5 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_5 = box_reduced_single_err
                        ad[0].BPM_FIBER_5 = bpm_single_fiber
                else:
                    """
                    Dealing with the case that we have no optimal extraction, so the reduced order
                    is the box extraction as opposed to the optimal extraction.
                    """
                    box_reduced_single_fiber_order_key = np.array(list(
                        box_reduced_stripes[f].keys()), dtype=float)
                    box_reduced_single_fiber = np.array(list(
                        box_reduced_stripes[f].values()), dtype=float)
                    box_reduced_single_err = np.array(list(
                        box_reduced_err[f].values()), dtype=float)
                    bpm_single_fiber = np.array(list(
                        extracted_bpms[f].values()), dtype=int)
                    # Update the extensions based on which fiber we have.  Again, probably a cleaner way to do this.
                    if f == 'fiber_1':
                        ad[0].REDUCED_ORDERS_FIBER_1 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_1 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_1 = box_reduced_single_err
                        ad[0].BPM_FIBER_1 = bpm_single_fiber
                    if f == 'fiber_2':
                        ad[0].REDUCED_ORDERS_FIBER_2 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_2 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_2 = box_reduced_single_err
                        ad[0].BPM_FIBER_2 = bpm_single_fiber
                    if f == 'fiber_3':
                        ad[0].REDUCED_ORDERS_FIBER_3 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_3 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_3 = box_reduced_single_err
                        ad[0].BPM_FIBER_3 = bpm_single_fiber
                    if f == 'fiber_4':
                        ad[0].REDUCED_ORDERS_FIBER_4 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_4 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_4 = box_reduced_single_err
                        ad[0].BPM_FIBER_4 = bpm_single_fiber
                    if f == 'fiber_5':
                        ad[0].REDUCED_ORDERS_FIBER_5 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_5 = box_reduced_single_fiber
                        ad[0].BOX_REDUCED_ERR_5 = box_reduced_single_err
                        ad[0].BPM_FIBER_5 = bpm_single_fiber

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
        Utilized in the dynamic wavelength calibration recipe as we do 
        not have any flats available at this stage. 

        Parameters
        ----------
        adinputs with STRIPES, F_STRIPES, and STRIPES_MASKS 'extensions' as
            dicts of sparse arrays
        
        Returns
        -------
        adinputs with box extracted orders for each fiber as
        well as uncertainties calculated during the box extraction
        """
        
    @staticmethod
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

        """
        Create matrix containing all indices for a given slit height.
        For the y matrix, the matrix is values from - slit height to + slit height,
        repeating nx times.  For the x matrix, it is values 0 to nx repeated 2* slit height.
        Both have dimensions (2*slit height) x nx.
        """
        
        slit_indices_y = np.arange(-slit_height, slit_height
                                   ).repeat(nx).reshape((2 * slit_height, nx))
        slit_indices_x = np.tile(np.arange(nx), 2 * slit_height
                                 ).reshape((2 * slit_height, nx))

        indices = np.rint(slit_indices_y + y).astype(int)
        valid_indices = np.logical_and(indices < ny, indices > 0)

        """
        Create sparse matrix of dimensions ny x nx, with row indices given by indices[valid_indices] 
        and column indices sit_indices_x[valid_indices].
        """
        mat = sparse.coo_matrix((data[indices[valid_indices],
            slit_indices_x[valid_indices]],(indices[valid_indices],
            slit_indices_x[valid_indices])), shape=(ny, nx))
        return mat.tocsc()

    def _optimal_extraction_single_stripe(self, stripe=None, flat_stripe=None,
                                          gain=2.7, read_noise=1.23,
                                          back_var=None, mask=None,
                                          full_output=False,
                                          s_clip=5.0, penalty=1.0):
        """
        Performs optimal extraction of a single stripe.
        Based on the algorithm described by Horne et al. 1986, PASP, 98, 609.

        Parameters
        ----------
            stripe (scipy.sparse.spmatrix): sparse matrix stripe from science
                frame to be extracted
            flat_stripe (scipy.sparse.spmatrix): sparse matrix stripe from
                flat frame to be used as profile
            gain (float): detector gain factor (conversion photons -> DN,
                given in e-/DN )
            read_noise (float): typical detector read noise given as the
                variance in DN
            back_var (np.ndarray): background variance (from scattered light
                model or dark, otherwise 0)
            mask (np.ndarray): bad pixel mask. Note: the mask will be modified
            by iterative algorithm that looks for outliers during optimal
                extraction
            full_output (bool): if True, returns all intermediate
                results for debugging/testing
            s_clip (float): sigma clipping value for optimal extraction
            penalty (float): scaling factor for global per-order profile
                mismatch correction between expected (flat) and found (science).
                Set 0 for no correction

        Returns
        -------
            tuple(np.ndarray, np.ndarray, dict): (optimal extracted spectrum,
            box extracted spectrum, dict of additional intermediate results
            if full_output was True)

        """
        log = self.log
        # fix back_var
        if mask is None:
            mask = stripe.copy()
            mask[:, :] = 1
            back_var = stripe.copy()
            back_var[:, :] = 0
        else:
            # copy mask, because it will be modified
            mask = mask.copy()
            back_var = stripe.copy()
            back_var[:, :] = 0  # back_var.copy()

        # box extracted spectrum
        stand_spec0 = self._box_extract_single_stripe(stripe, mask)
        stand_spec_err = np.sqrt(stand_spec0/gain)
        stand_spec = stand_spec0.copy()

        # flat data
        box_extracted_flat_stripe = self._box_extract_single_stripe(
            flat_stripe, mask)  # spatial sums for flat along disp

        # cut stripe sparse matrix into numpy array
        # find the spatial columns utilized along entire stripe
        # (greater than slit height because of stripe path)
        sparse_vcols = np.array(~np.all(stripe.todense() == 0, axis=1)
                                ).reshape(-1)
        # strip stripe to the inclusive nonzero rows
        stripe = np.array(stripe.todense()[sparse_vcols])
        # strip other sparse matricies similarily
        mask = np.array(mask.todense()[sparse_vcols])
        back_var = np.array(back_var.todense()[sparse_vcols])
        flat_stripe = np.array(flat_stripe.todense()[sparse_vcols])
        # use ~middle column slit height as slit height pass
        sparse_vrows = np.count_nonzero((stripe != 0).T[1500])
        data = np.zeros((sparse_vrows, stripe.shape[1]))
        diff_save = data.copy()
        new_mask = data.copy()  # create actual limit numpy arrays
        new_back_var = data.copy()  # create actual limit numpy arrays
        profile = data.copy()
        # for column in stripe
        for i in np.arange(stripe.shape[1]):
            # if column is slit height (not edge of chip)
            if np.nonzero((stripe != 0).T[i])[0].shape[0] == data.shape[0]:
                # if the flat stripe aka profile is defined in that column
                if box_extracted_flat_stripe[i] > 1E-12:
                    # limit the column to the nonzero data points
                    data[:, i] = stripe[np.nonzero((stripe != 0).T[i])[0], i]
                    new_mask[:, i] = mask[np.nonzero((stripe != 0).T[i])[0], i]
                    new_back_var[:, i] = back_var[np.nonzero(
                        (stripe != 0).T[i])[0], i]
                    profile[:, i] = flat_stripe[np.nonzero(
                        (stripe != 0).T[i])[0], i] / box_extracted_flat_stripe[
                        i]
            else:
                # cancel iteration on column
                new_mask[:, i] = 0
                new_back_var[:, i] = 0

        # if sum is less than zero, whole column is cancelled for this stripe
        new_mask[:, stand_spec0 < 1E-12] = 0
        # update original variables
        mask = new_mask.copy()
        back_var = new_back_var.copy() * mask
        stripe = data
        data_var = abs(stripe.copy()) / gain + back_var + read_noise

        # final output
        flux = np.zeros(len(stand_spec0))
        var = flux.copy()

        # Calculate a first guess of the difference between stripe and profile
        expected = profile * mask * stand_spec
        actual = stripe * mask
        diff = actual - expected

        # Calculate the median of the difference along the order.
        # Represents the 'global' mismatch between flat and science profiles
        # and helps correct for the 'drift' problem in x-dispersion.
        diff_aver = penalty * np.abs(median_filter(diff, size=(1, 201)))

        # avoid already caught bad pixels
        # aka where whole column is zero in bpm or stripe is too small
        good_disp = np.nonzero(np.array(~np.any(mask == 0, axis=0)))[0]
        reject_tracker = np.zeros_like(stand_spec0)

        for h in good_disp:
            # in good columns
            # flat column with mask scaled to data total
            expected = profile[:, h] * mask[:, h] * stand_spec[h]
            # actual column data with mask
            actual = stripe[:, h] * mask[:, h]
            diff = actual - expected

            data_var[:, h] = abs(stand_spec[h] * profile[:, h]) / gain + \
                             back_var[:, h] + read_noise
            noise_rev = 1 / np.sqrt(data_var[:, h])
            diff_save[:, h] = diff * noise_rev
            # track rejected mismatch between profiles based on s_clip value
            reject_index = np.nonzero(((np.abs(diff) - np.abs(diff_aver[:, h]))
                                       * noise_rev) >= s_clip)[0]

            while len(reject_index) > 0:
                # find and remove the worst mismatch
                worst = np.argmax((np.abs(diff) - np.abs(diff_aver[:, h]))
                                  * noise_rev)
                reject_tracker[h] = reject_tracker[h] + 1
                mask[worst, h] = 0
                # (re-)calculate the optimal extraction
                denom = np.sum(profile[:, h] * profile[:, h] * mask[:, h]
                               / data_var[:, h])
                stand_spec[h] = np.sum(profile[:, h] * mask[:, h] * stripe[:, h]
                                       / data_var[:, h]) / denom

                expected = profile[:, h] * mask[:, h] * stand_spec[h]
                actual = stripe[:, h] * mask[:, h]
                diff = actual - expected

                data_var[:, h] = abs(stand_spec[h] * profile[:, h]) / gain + \
                                 back_var[:, h] + read_noise
                noise_rev = 1 / np.sqrt(data_var[:, h])
                # diff_save[:,h] = diff*noise_rev
                reject_index = np.nonzero(((np.abs(diff) - np.abs(
                    diff_aver[:, h])) * noise_rev) >= s_clip)[0]
                # if more than half the pixels in a column have been rejected,
                # reject hte whole column
                if np.count_nonzero(actual[3:-5]) < len(profile[3:-5, h]) / 2.:
                    reject_index = np.array([])
                    mask[:, h] = 0
                    flux[h] = 0
                    log.warning(f'Too many bad pixels in column {h}, reject column')
            # if the column is still considered good after all the rejections
            # are complete, calculate a final optimal extraction
            if np.count_nonzero(mask[:, h]) > 0:
                denom = np.sum(profile[:, h] * profile[:, h] * mask[:, h] /
                               data_var[:, h])
                flux[h] = np.sum(profile[:, h] * mask[:, h] * stripe[:, h] /
                                 data_var[:, h]) / denom
                var[h] = np.sum(profile[:, h] * mask[:, h]) / denom
        # from the optimal extraction, values of 0 in the flux are bad pixels
        flux[flux == 0] = np.nan
        var[var == 0] = np.nan
        stand_spec0[stand_spec0 == 0] = np.nan
        # check to see if the total rejected pixels in a stripe is significant
        total_count = np.sum(reject_tracker > 0)
        # rejection fraction limit, found by instrument team
        substantial_count = np.sum(
            np.abs(stand_spec0 - stand_spec)[reject_tracker > 0] > 0.005 *
            (stand_spec0[reject_tracker > 0]))
        irrelevant_count = total_count - substantial_count

        if total_count > 0:
            log.fullinfo(
                f'Rejected {np.sum(reject_tracker):.0f} pixels in ' +
                f'{total_count} columns during optimal extraction.')
            irrelevant_count_percentage = irrelevant_count / total_count * 100
            if irrelevant_count_percentage > 30:
                log.warning(
                    f'Rejections with flux changes < 0.5%: {irrelevant_count}' +
                    f'({irrelevant_count_percentage:.0f}%)')
            else:
                log.info(
                    f'Rejections with flux changes < 0.5%: {irrelevant_count}' +
                    f'({irrelevant_count_percentage:.0f}%)')
        else:
            log.info('No rejections')

        # return all intermediate results for debugging/testing
        if full_output:
            return flux, var, stand_spec0, stand_spec_err, \
                {'noise': var, 'acceptancemask': mask, "stripe": stripe,
                                            "initial_sigma": diff_save}
        return flux, var, stand_spec0, stand_spec_err, {'noise': var}

    @staticmethod
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
