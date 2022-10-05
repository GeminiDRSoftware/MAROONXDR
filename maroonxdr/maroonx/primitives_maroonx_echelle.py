#
#                                                                       DRAGONS
#
#                                                 primitives_maroonx_echelle.py
# ------------------------------------------------------------------------------

from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.gemini import gemini_tools as gt
import numpy as np
from scipy.ndimage import median_filter
from scipy import ndimage
import copy
import pandas as pd
from astropy.table import Table, vstack
import scipy.sparse as sparse
from geminidr.core import Spect
from .primitives_maroonx import MAROONX
from . import parameters_maroonx_echelle
import matplotlib.pyplot as plt
from recipe_system.utils.decorators import parameter_override
import astrodata  # replace with correct primative call for caldb get flat & dark
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

    def standardizeWCS(self, adinputs=None, suffix=None, **params):
        """
        MAROONXDR-specific version of validateData to ignore the invalid WCS
        exception.
        """
        try:
            super().standardizeWCS(adinputs, suffix=suffix)
        except TypeError as e:
            if 'The value must be' not in str(e):
                raise
        return adinputs

    def extractStripes(self, adinputs=None, flat=None, dark=None, skip_dark=None, slit_height=10, debug_level=0, **params):
        """
        Extracts the stripes from the original 2D spectrum to a sparse array, containing only relevant pixels.

        This function marks all relevant pixels for extraction. Using the provided dictionary P_id it iterates over all
        stripes in the image and saves a sparse matrix for each stripe.

        Args:
            img (np.ndarray): 2d echelle spectrum
            p_id (Union[dict, str]): dictionary as returned by :func:`~identify_stripes` or path to file
            slit_height (int): total slit height in px
            debug_level (int): debug level

        Returns:
            dict: dictionary of the form {fiber_number:{order: scipy.sparse_matrix}}

        """
        if skip_dark is None:
            skip_dark = [1, 2, 3, 4, 5]
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            if flat is None:  # fix so self.caldb.get_processed_flat(ad) call works
                if 'BLUE' in ad.tags:
                    flat_ad = astrodata.open('calibrations/processed_flat/20200911T220106Z_FDDDF_b_0002_FFFFF_flat.fits')
                elif 'RED' in ad.tags:
                    flat_ad = astrodata.open('calibrations/processed_flat/20200911T220106Z_FDDDF_r_0000_FFFFF_flat.fits')
                # flat_ad = self.caldb.get_calibrations(ad,caltype='processed_flat')
                # flat_ad = astrodata.open(flat)
                else:
                    log.warning("No extraction will be made on {}, since no "
                                "flat was found or specified".format(ad.filename))

            if dark is None:
                if 'BLUE' in ad.tags:
                    dark_ad = astrodata.open('calibrations/processed_dark/20201120T165414Z_DDDDE_b_0900_dark.fits')
                elif 'RED' in ad.tags:
                    dark_ad = astrodata.open('calibrations/processed_dark/20201120T165414Z_DDDDE_r_0900_dark.fits')
                # dark_ad = self.caldb.get_calibrations(ad,caltype='processed_dark')
                # dark_ad = astrodata.open(flat)
                else:
                    log.warning("No dark subtraction will be made to {} prior to stripe extraction, since no "
                                "dark was found or specified".format(ad.filename))

            stripes = {}
            f_stripes = {}
            p_id = flat_ad[0].STRIPES_ID
            if len(p_id) != 30:
                log.warning("No extraction will be made on {}, since "
                            "flat in use is incorrectly formated".format(ad.filename))

            ad[0].STRIPES_ID = p_id  # record used stripe reference to science frame
            #repackage p_id as dict
            p_id = {'fiber_1': dict((colname, p_id[0:6][colname].data) for colname in p_id.colnames),
                        'fiber_2': dict((colname, p_id[6:12][colname].data) for colname in p_id.colnames),
                        'fiber_3': dict((colname, p_id[12:18][colname].data) for colname in p_id.colnames),
                        'fiber_4': dict((colname, p_id[18:24][colname].data) for colname in p_id.colnames),
                        'fiber_5': dict((colname, p_id[24:][colname].data) for colname in p_id.colnames)}

            for f in p_id.keys():  # read stripes for both science frame and flat frame
                adint = copy.deepcopy(ad)
                flatint = copy.deepcopy(flat_ad)
                if int(f[-1]) in skip_dark:
                    log.fullinfo(f'No dark subtracted for fiber {f[-1]}')
                    adint.data[0] = ad.data[0]
                else:
                    adint.data[0] = ad.data[0] - dark_ad.data[0]
                for o, p in p_id[f].items():
                    stripe = self._extract_single_stripe(adint, p, slit_height, debug_level)
                    f_stripe = self._extract_single_stripe(flatint, p, slit_height, debug_level)
                    if f in stripes:
                        stripes[f].update({o: stripe})
                        f_stripes[f].update({o: f_stripe})
                    else:
                        stripes[f] = {o: stripe}
                        f_stripes[f] = {o: f_stripe}

            ad[0].STRIPES = stripes
            ad[0].F_STRIPES = f_stripes  # could find a way to make the flat info retrieve more efficient
            ad.update_filename(suffix=params['suffix'], strip=True)
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def optimalExtraction(self, adinputs=None, opt_extraction=None, back_var=None, full_output=False,
                           penalty=None, s_clip=None, debug_level=0, **params):
        """
        Optimal extraction of the 2d echelle spectrum.

        This function performs an optimal extraction of a 2d echelle spectrum.
        A given flat field spectrum is used to generate normalized 'profiles' that are used as weighting functions for the
        spectrum that is going to be extracted.
        The algorithm further checks for outliers and rejects them. This is to prevent contributions from cosmic hits.

        Args:
            stripes and flat stripes found in adinputs
            stripes (dict): stripes to be extracted. As returned by :func:`~extract_stripes`
            flat_stripes (dict): corresponding stripes of a (master) flat field
            opt_extraction (list): list with fibers considered for optimal extraction
            gain (float): CCD gain, gotten from astrodata object
            read_noise (float): CCD read noise, gotten from astrodata object
            back_var ():
            mask (np.ndarray): binary pixel mask. Only pixel where mask==1 will be considered for extraction
            full_output (bool): if True, an additional set of intermediate products will be returned / saved
            debug_level (int): debug level

        Returns:
            tuple(dict, dict, dict): dictionary of optimal and box reduced stripes and optimal extraction intermediate
            results for debugging
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        if opt_extraction is None:
            opt_extraction = [0,2,3,4,0]
        optimal_reduced_stripes = {}
        box_reduced_stripes = {}
        optimal_reduced_err = {}
        optimal_reduced_2d_arrays = {}
        for ad in adinputs:
            ad[0].REDUCED_ORDERS_FIBER_1 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_FIBER_1 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_ERR_1 = np.zeros(shape=[1,1])
            ad[0].BOX_REDUCED_FIBER_1 = np.zeros(shape=[1,1])
            ad[0].REDUCED_ORDERS_FIBER_2 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_FIBER_2 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_ERR_2 = np.zeros(shape=[1,1])
            ad[0].BOX_REDUCED_FIBER_2 = np.zeros(shape=[1,1])
            ad[0].REDUCED_ORDERS_FIBER_3 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_FIBER_3 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_ERR_3 = np.zeros(shape=[1,1])
            ad[0].BOX_REDUCED_FIBER_3 = np.zeros(shape=[1,1])
            ad[0].REDUCED_ORDERS_FIBER_4 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_FIBER_4 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_ERR_4 = np.zeros(shape=[1,1])
            ad[0].BOX_REDUCED_FIBER_4 = np.zeros(shape=[1,1])
            ad[0].REDUCED_ORDERS_FIBER_5 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_FIBER_5 = np.zeros(shape=[1,1])
            ad[0].OPTIMAL_REDUCED_ERR_5 = np.zeros(shape=[1,1])
            ad[0].BOX_REDUCED_FIBER_5 = np.zeros(shape=[1,1])
            stripes = ad[0].STRIPES
            gain = ad.gain()[0][0] # fix
            read_noise = ad.read_noise()[0][0] # fix
            mask = None  # get bpm
            flat_stripes = ad[0].F_STRIPES
            for f in stripes.keys():
                if int(f[-1]) in opt_extraction:
                    for o, stripe in stripes[f].items():
                        log.fullinfo(f'Optimum extraction in {f}, order {o}')
                        flux, var, stand_spec, fo = self._optimal_extraction_single_stripe(stripe, flat_stripes[f][o],
                                                                                     gain=gain,
                                                                                     read_noise=read_noise,
                                                                                     back_var=back_var,
                                                                                     mask=mask,
                                                                                     s_clip=s_clip,
                                                                                     penalty=penalty,
                                                                                     debug_level=debug_level,
                                                                                     full_output=full_output)

                        if f in optimal_reduced_stripes:
                            optimal_reduced_stripes[f].update({o: flux})
                            optimal_reduced_err[f].update({o: var})
                            optimal_reduced_2d_arrays[f].update({o: fo})
                            box_reduced_stripes[f].update({o: stand_spec})
                        else:
                            optimal_reduced_stripes[f] = {o: flux}
                            optimal_reduced_err[f] = {o: var}
                            optimal_reduced_2d_arrays[f] = {o: fo}
                            box_reduced_stripes[f] = {o: stand_spec}
                else:
                    for o, stripe in stripes[f].items():
                        log.fullinfo(f'Only box extraction in {f}, order {o}')
                        stand_spec = self._box_extract_single_stripe(stripe, mask)

                        if f in box_reduced_stripes:
                            box_reduced_stripes[f].update({o: stand_spec})
                        else:
                            box_reduced_stripes[f] = {o: stand_spec}
            for f in stripes.keys():
                if f in optimal_reduced_stripes.keys():
                    optimal_reduced_single_fiber = np.array(list(optimal_reduced_stripes[f].values()),dtype=int)
                    optimal_reduced_single_fiber_order_key = np.array(list(optimal_reduced_stripes[f].keys()),dtype=float)
                    optimal_reduced_single_fiber_err = np.array(list(optimal_reduced_err[f].values()),dtype=float)
                    box_reduced_stripes_single_fiber = np.array(list(box_reduced_stripes[f].values()),dtype=float)
                    if f == 'fiber_1':
                        ad[0].REDUCED_ORDERS_FIBER_1 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_1 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_1 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_1 = box_reduced_stripes_single_fiber
                    if f == 'fiber_2':
                        ad[0].REDUCED_ORDERS_FIBER_2 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_2 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_2 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_2 = box_reduced_stripes_single_fiber
                    if f == 'fiber_3':
                        ad[0].REDUCED_ORDERS_FIBER_3 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_3 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_3 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_3 = box_reduced_stripes_single_fiber
                    if f == 'fiber_4':
                        ad[0].REDUCED_ORDERS_FIBER_4 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_4 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_4 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_4 = box_reduced_stripes_single_fiber
                    if f == 'fiber_5':
                        ad[0].REDUCED_ORDERS_FIBER_5 = optimal_reduced_single_fiber_order_key
                        ad[0].OPTIMAL_REDUCED_FIBER_5 = optimal_reduced_single_fiber
                        ad[0].OPTIMAL_REDUCED_ERR_5 = optimal_reduced_single_fiber_err
                        ad[0].BOX_REDUCED_FIBER_5 = box_reduced_stripes_single_fiber
                else:
                    box_reduced_single_fiber_order_key = np.array(list(box_reduced_stripes[f].keys()),dtype=float)
                    box_reduced_stripes_single_fiber = np.array(list(box_reduced_stripes[f].values()), dtype=float)
                    if f == 'fiber_1':
                        ad[0].REDUCED_ORDERS_FIBER_1 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_1 = box_reduced_stripes_single_fiber
                    if f == 'fiber_2':
                        ad[0].REDUCED_ORDERS_FIBER_2 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_2 = box_reduced_stripes_single_fiber
                    if f == 'fiber_3':
                        ad[0].REDUCED_ORDERS_FIBER_3 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_3 = box_reduced_stripes_single_fiber
                    if f == 'fiber_4':
                        ad[0].REDUCED_ORDERS_FIBER_4 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_4 = box_reduced_stripes_single_fiber
                    if f == 'fiber_5':
                        ad[0].REDUCED_ORDERS_FIBER_5 = box_reduced_single_fiber_order_key
                        ad[0].BOX_REDUCED_FIBER_5 = box_reduced_stripes_single_fiber
            del ad[0].STRIPES  # do not want to save the sparse matricies
            del ad[0].F_STRIPES
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=False)
        return adinputs

    @staticmethod
    def _extract_single_stripe(ad=None, polynomials=None, slit_height=10, debug_level=0):
        """
        Extracts single stripe from 2d image.

        This function returns a sparse matrix containing all relevant pixel for a single stripe for a given polynomial p
        and a given slit height.

        Args:
            polynomials (np.ndarray): polynomial coefficients
            img (np.ndarray): 2d echelle spectrum
            slit_height (int): total slit height in pixel to be extracted
            debug_level (int): debug level

        Returns:
            scipy.sparse.csc_matrix: extracted spectrum

        """
        ny, nx = ad.data[0].shape
        if isinstance(slit_height, np.ndarray):
            slit_height = int(slit_height[0])

        xx = np.arange(nx)
        y = np.poly1d(polynomials)(xx)

        slit_indices_y = np.arange(-slit_height, slit_height).repeat(nx).reshape((2 * slit_height, nx))
        slit_indices_x = np.tile(np.arange(nx), 2 * slit_height).reshape((2 * slit_height, nx))

        indices = np.rint(slit_indices_y + y).astype(int)
        valid_indices = np.logical_and(indices < ny, indices > 0)

        if debug_level > 3:
            plt.figure()
            plt.imshow(ad.data[0], origin='lower')
            ind_img = np.zeros_like(ad.data[0])
            ind_img[indices[valid_indices], slit_indices_x[valid_indices]] = 1
            plt.imshow(ind_img, alpha=0.5, origin='lower')
            plt.show()

        mat = sparse.coo_matrix((ad.data[0][indices[valid_indices], slit_indices_x[valid_indices]],
                                 (indices[valid_indices], slit_indices_x[valid_indices])), shape=(ny, nx))
        return mat.tocsc()

    def _optimal_extraction_single_stripe(self, stripe=None, flat_stripe=None, gain=2.7, read_noise=1.23, back_var=None, mask=None,
                                         debug_level=0, full_output=False, s_clip=5.0, penalty=1.0):
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

        Returns:
            tuple(np.ndarray, np.ndarray, dict): (optimal extracted spectrum, box extracted spectrum, dict of additional
            intermediate results if full_output was True)

        """
        log = self.log
        # back_var_constant = back_var
        if mask is None:
            mask = stripe.copy()
            mask[:, :] = 1
            back_var = stripe.copy()
            back_var[:, :] = 0
        else:
            # copy mask, because it will be modified
            mask = mask.copy()
            back_var = back_var.copy()

        # box extracted spectrum
        stand_spec0 = self._box_extract_single_stripe(stripe, mask)  # direct box extraction,
        stand_spec = stand_spec0.copy()

        # flat data
        box_extracted_flat_stripe = self._box_extract_single_stripe(flat_stripe, mask)  # spatial sums for flat along disp

        # cut stripe sparse matrix into numpy array
        # find the spatial columns utilized along entire stripe (greater than slit height because of stripe path)
        sparse_vcols = np.array(~np.all(stripe.todense() == 0, axis=1)).reshape(-1)
        stripe = np.array(stripe.todense()[sparse_vcols])  # strip stripe to the inclusive nonzero rows
        mask = np.array(mask.todense()[sparse_vcols])  # strip mask similarly
        back_var = np.array(back_var.todense()[sparse_vcols])  # strip background variance map similarly
        flat_stripe = np.array(flat_stripe.todense()[sparse_vcols])
        sparse_vrows = np.count_nonzero((stripe != 0).T[1500])  # use ~middle column slit height as slit height pass
        data = np.zeros((sparse_vrows, stripe.shape[1]))
        diff_save = data.copy()
        new_mask = data.copy()  # create actual limit numpy arrays
        new_back_var = data.copy()  # create actual limit numpy arrays
        profile = data.copy()
        for i in np.arange(stripe.shape[1]):
            if np.nonzero((stripe != 0).T[i])[0].shape[0] == data.shape[0]:  # if column is slit height (not edge of chip)
                if box_extracted_flat_stripe[i] > 1E-12:
                    data[:, i] = stripe[np.nonzero((stripe != 0).T[i])[0], i]  # write data
                    new_mask[:, i] = mask[np.nonzero((stripe != 0).T[i])[0], i]
                    new_back_var[:, i] = back_var[np.nonzero((stripe != 0).T[i])[0], i]
                    profile[:, i] = flat_stripe[np.nonzero((stripe != 0).T[i])[0], i] / box_extracted_flat_stripe[i]
            else:
                new_mask[:, i] = 0  # could be optimized
                new_back_var[:, i] = 0  # could be optimized

        new_mask[:, stand_spec0 < 1E-12] = 0  # if sum is less than zero, whole column is cancelled for this stripe
        mask = new_mask.copy()
        back_var = new_back_var.copy() * mask
        stripe = data  # return variables 'mask' and 'stripe' to the naming conventions
        data_var = abs(stripe.copy()) / gain + back_var + read_noise

        # final output
        flux = np.zeros(len(stand_spec0))
        var = flux.copy()

        # Calculate a first guess of the difference between stripe and scaled flat_stripe
        expected = profile * mask * stand_spec
        actual = stripe * mask
        diff = actual - expected

        # Calculate the median of the difference along the order. This represents the 'global' mismatch between flat and science profile
        # and helps correct for the 'drift' problem in x-dispersion.
        diff_aver = penalty * np.abs(median_filter(diff, size=(1, 201)))

        # # avoid already caught bad pixels (whole column is zero in bpm or stripe is too small)
        good_disp = np.nonzero(np.array(~np.any(mask == 0, axis=0)))[0]
        reject_tracker = np.zeros_like(stand_spec0)

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
                if debug_level >= 3:
                    log.fullinfo(f'Outlier found in column {h}, pixel {worst}')
                    fig, ax = plt.subplots(2, 1)
                    actual_plot = actual.copy()
                    expected_plot = expected.copy()
                    actual_plot[actual_plot == 0] = np.nan
                    expected_plot[expected_plot == 0] = np.nan
                    ax[0].plot(actual_plot, 'r')
                    ax[0].plot(expected_plot, 'b')
                    ax[1].plot(np.abs(diff), 'r')
                    ax[1].plot(np.abs(diff) - np.abs(diff_aver[:, h]), 'g')
                    ax[1].plot(np.sqrt(data_var[:, h]) * s_clip, 'b')
                    plt.show()

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

                if np.count_nonzero(actual[3:-5]) < len(profile[3:-5, h]) / 2.:
                    reject_index = np.array([])
                    mask[:, h] = 0
                    flux[h] = 0
                    log.warning(f'Too many bad pixels in column {h}, reject column')
            if np.count_nonzero(mask[:, h]) > 0:
                denom = np.sum(profile[:, h] * profile[:, h] * mask[:, h] / data_var[:, h])
                flux[h] = np.sum(profile[:, h] * mask[:, h] * stripe[:, h] / data_var[:, h]) / denom
                var[h] = np.sum(profile[:, h] * mask[:, h]) / denom

        flux[flux == 0] = np.nan
        var[var == 0] = np.nan

        stand_spec0[stand_spec0 == 0] = np.nan

        total_count = np.sum(reject_tracker > 0)
        substantial_count = np.sum(
            np.abs(stand_spec0 - stand_spec)[reject_tracker > 0] > 0.005 * (stand_spec0[reject_tracker > 0]))
        irrelevant_count = total_count - substantial_count

        if total_count > 0:
            log.fullinfo(
                f'Rejected {np.sum(reject_tracker):.0f} pixels in {total_count} columns during optimal extraction.')
            irrelevant_count_percentage = irrelevant_count / total_count * 100
            if irrelevant_count_percentage > 30:
                log.warning(
                    f'Rejections with flux changes < 0.5%: {irrelevant_count} ({irrelevant_count_percentage:.0f}%)')
            else:
                log.info(
                    f'Rejections with flux changes < 0.5%: {irrelevant_count} ({irrelevant_count_percentage:.0f}%)')
        else:
            log.info('No rejections')

        if debug_level >= 2:
            # fig, ax = plt.subplots(3, 1)
            fig, ax = plt.subplots(3, 1, sharex='all')
            ax[0].imshow(np.vstack((stripe, np.zeros_like(stripe), profile * flux)), interpolation='none')
            ax[0].scatter(np.where(mask == 0)[1], np.where(mask == 0)[0], c='r', marker='+')

            ax[1].plot(flux - stand_spec0, 'r', label='Difference box vs optimal')
            ax[2].plot(stand_spec0, 'b', label='Box extraction')
            ax[2].plot(flux, 'g', label='Optimal extraction')
            plt.legend()

            plt.show()

        # return all intermediate results for debugging/testing
        if full_output:
            return flux, var, stand_spec0, {'noise': var, 'acceptancemask': mask, "stripe": stripe,
                                            "initial_sigma": diff_save}  # {'noise': var, 'profile': profile, 'rejectionmask': sparse.csr_matrix(mask),
            # 'expected': expected, 'actual': actual}
        else:
            return flux, var, stand_spec0, {'noise': var}

    @staticmethod
    def _box_extract_single_stripe(stripe=None, mask=None):
        """
        Box extraction of a single stripe.

        Args:
            stripe (sparse.matrix): stripe to be extracted
            mask (np.ndarray): binary pixel mask. Where mask==1: values will contribute. Same shape as stripe.toarray()

        Returns:
            np.ndarray: box extracted flux
        """
        if mask is None:
            mask = stripe.copy()
            mask.data[:] = 1.
        return np.array(np.sum(stripe.multiply(mask), axis=0).T).flatten()
