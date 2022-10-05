#
#                                                                       DRAGONS
#
#                                                         primitives_maroonx.py
# ------------------------------------------------------------------------------
import os

import pandas as pd
from gempy.gemini import gemini_tools as gt
import numpy as np
import copy
from scipy import ndimage
from scipy.ndimage import median_filter
from astropy.table import Table, vstack
import scipy.sparse as sparse
from geminidr.gemini.primitives_gemini import Gemini
from geminidr.core import CCD, NearIR, primitives_preprocess
from astropy.io import fits
from astropy.stats import SigmaClip
from . import parameters_maroonx
from photutils import Background2D, MedianBackground
from .lookups import timestamp_keywords as maroonx_stamps
from .lookups import siddb as maroonx_siddb
from .lookups import maskdb as maroonx_maskdb
#from gempy.library import astromodels
from recipe_system.utils.decorators import parameter_override


import matplotlib.pyplot as plt
# ------------------------------------------------------------------------------

@parameter_override
class MAROONX(Gemini, CCD, NearIR):
    """
    This class inherits from the level above.  Any primitives specific
    to MAROON-X can go here. 
    """

    tagset = {"GEMINI", "MAROONX"}

    def __init__(self, adinputs, **kwargs):
        super(MAROONX, self).__init__(adinputs, **kwargs)
        self._param_update(parameters_maroonx)
        # Add MAROON-X specific timestamp keywords
        self.timestamp_keys.update(maroonx_stamps.timestamp_keys)

    def checkArm(self, adinputs=None, **params):
        """
        Check that the the camera arm is consistent through all input files, first file sets

        Parameters
        ----------
        adinputs

        Returns
        -------
        adinputs that pass test

        """
        log = self.log
        arm_set = 'BLUE' if 'BLUE' in adinputs[0].tags else 'RED' if 'RED' in adinputs[0].tags else 'UNDEFINED'
        adoutputs = []
        for ad in adinputs:
            if arm_set not in ad.tags:
                log.warning("Not all frames taken with the same camera arm, restricting set to first arm used in list")
                log.warning('Not analyzing frame: '+ad.filename)
            else:
                ad.update_filename(suffix=params['suffix'], strip=True)
                adoutputs.append(ad)
        return adoutputs

    def validateData(self, adinputs=None, suffix=None):
        """
        MAROONXDR-specific version of validateData to ignore the invalid WCS
        exception.
        """
        try:
            super().validateData(adinputs, suffix=suffix)
        except ValueError as e:
            if 'valid WCS' not in str(e):
                raise
        return adinputs
    def correctImageOrientation(self, adinputs=None, debug_level=0, **params):
        """
        Correct image orientation to proper echelle format.

        flips image so that left lower corner is bluest wavelength, upper right corner is reddest wavelength.
        Echelle orders go from left to right.

        Args:
            adinputs --- previously --- img (np.ndarray): input 2d echelle spectrum
            debug_level (int): debug level
        Returns:
            np.ndarray: image with correct orientation
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        adoutputs = []
        for ad in adinputs:
            adout = copy.deepcopy(ad)
            if debug_level > 0:
                plt.figure()
                plt.title('Original Image')
                plt.imshow(ad.data[0], origin='lower')
                plt.show()

            if ad.image_orientation()['vertical orientation flip'] and ad.image_orientation()['horizontal orientation flip']:  # flip up-down
                adout[0].data = np.fliplr(np.flipud(ad[0].data))

            if debug_level > 0:
                plt.figure()
                plt.title('Orientation Corrected image')
                plt.imshow(adout.data[0], origin='lower')
                plt.show()
            adout.update_filename(suffix=params['suffix'], strip=True)
            adoutputs.append(adout)

        gt.mark_history(adoutputs, primname=self.myself(), keyword=timestamp_key)

        return adoutputs

    def checkND(self, adinputs=None, **params):
        """
        Check that the ND filter on the sim cal fiber is consistent through all input files, first file sets

        Parameters
        ----------
        adinputs

        Returns
        -------
        adinputs that pass test

        """
        log = self.log
        check_val = adinputs[0].filter_orientation()['ND']
        adoutputs = []
        if len(adinputs) > 1:
            for ad in adinputs:
                if check_val != ad.filter_orientation()['ND']:
                    log.warning("Not all frames have the same simcal ND filter setting, restricting set to first seen")
                else:
                    ad.update_filename(suffix=params['suffix'], strip=True)
                    adoutputs.append(ad)
            if len(adoutputs) == 1:
                raise IOError("Less than two frames found with first frame simcal ND filter setting")
        else:
            return adinputs
        return adoutputs

    def _get_sid_filename(self,ad):
        """
        Gets stripe id file for input frame

        Returns
        -------
        str/None: Filename of appropriate sid
        """
        log = self.log
        arm = ('b' if 'BLUE' in ad.tags else 'r')
        sid_dir = os.path.join(os.path.dirname(maroonx_siddb.__file__), 'SID')
        db_matches = sorted((k, v) for k, v in maroonx_siddb.sid_dict.items() if arm in k)
        if db_matches:
            sid = db_matches[-1][1]
        else:
            log.warning('No SID found for {}'.format(ad.filename))
            return None
        return sid if sid.startswith(os.path.sep) else os.path.join(sid_dir, sid)

    def findStripes(self,adinputs=None, deg_polynomial=5, median_filter=1, gauss_filter_sigma=3.5, min_peak=0.008,
                    debug_level=0,**params):
        """   Locates and fits stripes in a flat field or science echelle spectrum.

                Starting in the central column, the algorithm identifies peaks and traces each stripe to the edge of the detector
                by following the brightest pixels along each order. It then fits a polynomial to each stripe. To improve
                algorithm stability, the image is first median filtered and then smoothed with a gaussian. It not only eliminates
                noise, but also ensures that the cross section profile of the flat becomes peaked in the middle, which helps to
                identify the center of each stripe. Choose gauss_filter accordingly. To avoid false positives, only peaks above a
                certain (relative) intensity threshold are used.

                Args:
                    deg_polynomial (int): degree of the polynomial fit
                    median_filter (int): median filter parameter
                    gauss_filter_sigma (float): sigma of the gaussian filter used to smooth the image.
                    min_peak (float): minimum relative peak height
                    debug_level (int): debug level

        """
        # don't need to reference a lookup, old reference was just to default fit mins
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:

            ny, nx = ad.data[0].shape

            filtered_ad = ndimage.median_filter(ad.data[0], median_filter)
            filtered_ad = ndimage.gaussian_filter(filtered_ad, gauss_filter_sigma)

            # find peaks in center column
            data = filtered_ad[:, int(nx / 2)]
            peaks = np.r_[True, data[1:] >= data[:-1]] & np.r_[data[:-1] > data[1:], True]

            idx = np.logical_and(peaks, data > min_peak * np.max(data))
            maxima = np.arange(ny)[idx]

            # filter out maxima too close to the boundary to avoid problems
            maxima = maxima[maxima > 3]
            maxima = maxima[maxima < ny - 3]

            if debug_level > 2:
                plt.figure()
                plt.title('Local maxima')
                plt.plot(data)
                plt.plot(np.arange(ny)[peaks], data[peaks], 'r+')
                plt.plot(maxima, data[maxima], 'go')
                plt.show()

            n_order = len(maxima)
            log.fullinfo(f'Number of stripes found: {n_order}')

            orders = np.zeros((n_order, nx))
            # walk through to the left and right along the maximum of the order
            for m, row in enumerate(maxima):
                column = int(nx / 2)
                orders[m, column] = row
                start_row = row
                # walk right
                while column + 1 < nx:
                    column += 1
                    args = [start_row]
                    if start_row - 1 > 1:
                        args.append(start_row - 1)
                    else:
                        args.append(1)
                    if start_row + 1 < ny - 1:
                        args.append(start_row + 1)
                    else:
                        args.append(ny - 1)
                    # args = np.array(np.linspace(max(1, start_row - 1), min(start_row + 1, ny - 1), 3), dtype=int)
                    p = filtered_ad[args, column]
                    # new maximum
                    start_row = args[int(np.argmax(p))]
                    orders[m, column] = start_row

                column = int(nx / 2)
                start_row = row
                # walk left
                while column > 0:
                    column -= 1
                    args = [start_row]
                    if start_row - 1 > 1:
                        args.append(start_row - 1)
                    else:
                        args.append(1)
                    if start_row + 1 < ny - 1:
                        args.append(start_row + 1)
                    else:
                        args.append(ny - 1)

                    # args = np.array(np.linspace(max(1, start_row - 1), min(start_row + 1, ny - 1), 3), dtype=int)
                    # args = args[np.logical_and(args < ny, args > 0)]
                    p = filtered_ad[args, column]
                    # new maximum
                    start_row = args[int(np.argmax(p))]
                    orders[m, column] = start_row

            # do Polynomial fit for each order
            log.fullinfo(f'Fit polynomial of order {deg_polynomial} to each stripe')
            xx = np.arange(nx)
            polynomials = [np.poly1d(np.polyfit(xx, o, deg_polynomial)) for o in orders]
            if debug_level > 0:
                plt.figure()
                plt.title('Traced stripes')
                plt.imshow(filtered_ad, interpolation='none', vmin=np.min(ad.data[0]), vmax=0.5 * np.max(ad.data[0]),
                           cmap=plt.get_cmap('gray'), origin='lower')
                for p in polynomials:
                    plt.plot(xx, p(xx), 'g', alpha=0.5)
                plt.ylim((0, ny))
                plt.xlim((0, nx))
                plt.show()
            ad[0].STRIPES_LOC = polynomials
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=False)
        return adinputs

    def identifyStripes(self, adinputs=None, position_dir=None, selected_fibers=None, debug_level=0,**params):
        """
        Identifies the stripes and assigns the proper order and fiber number.

        Returns a dictionary {fiber_X: {order_N: polynomial_coefficients , order_M: polynomial_coefficients, ...}, ...}

        Args:
            adinputs (np.ndarray): 2D dark subtracted flat field image
            polynomials (list): list of polynomials as returned in flat at ad.nddata[0].meta['STRIPES_LOC']
            positions (array): lookup array of nominal y positions and fiber/order labels. Shape is Nx3, columns are
            [fibers, orders, y_positions], found in lookups/SID
            fibers (str): if not None,
            debug_level (int): debug level
            selected_fibers (list): fibers illuminated in the flat

        Returns:
            dict: dictionary with labeled (fiber/order) polynomial coefficients

        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        positions_dir = None
        for ad in adinputs:
            if position_dir is None:
                positions_dir = self._get_sid_filename(ad)
                log.info(positions_dir)
            else:
                positions_dir = position_dir
            with fits.open(positions_dir,'readonly') as f:
                positions = f[1].data
            ny, nx = ad.data[0].shape
            p_id = {}

            selected_fibers = list(np.asarray((selected_fibers.split(','))).astype(int))

            use_all_fibers = selected_fibers is None

            y_positions = positions['fiber_position'].astype('int')
            fibers = positions['identify_fiber'].astype('int')
            orders = positions['fiber_order'].astype('int')
            unique_fibers = np.unique(fibers)
            if selected_fibers is None:
                selected_fibers = unique_fibers

            for f in unique_fibers:
                if f not in selected_fibers:
                    idx = np.where(fibers == f)
                    y_positions = np.delete(y_positions, idx)
                    fibers = np.delete(fibers, idx)
                    orders = np.delete(orders, idx)

            used = np.zeros_like(y_positions)

            observed_y = []
            for i, p in enumerate(ad[0].STRIPES_LOC):
                observed_y.append(np.poly1d(p)(nx / 2))

            shifts = np.linspace(-200, 200, 2000)
            total_distances = []
            for shift in shifts:
                distances = []
                for y in observed_y:
                    closest_stripe_idx = np.argmin(np.abs(y_positions + shift - y))
                    distances.append(np.abs(y_positions[closest_stripe_idx] + shift - y))
                total_distances.append(np.array(distances).sum())

            # very important for correct fiber IDs on subsets of the illuminated fibers
            # weights towards minima closer to the initial guess
            if not use_all_fibers:
                total_distances += np.abs(shifts) * 2

            if debug_level > 2:
                plt.figure()
                plt.title("Total distance between observed stripe positions y and config stripe positions")
                plt.xlabel("y-shift [px]")
                plt.ylabel("Total distance")
                plt.plot(shifts, total_distances)
                plt.show()
            shift_calculated = shifts[np.argmin(total_distances)]

            for i, p in enumerate(ad[0].STRIPES_LOC):
                observed_y = np.poly1d(p)(nx / 2)
                closest_stripe_idx = np.argmin(np.abs(y_positions + shift_calculated - observed_y))
                if np.abs(y_positions[closest_stripe_idx] + shift_calculated - observed_y) < 7:
                    if used[closest_stripe_idx] == 0:
                        used[closest_stripe_idx] = 1
                        fiber = fibers[closest_stripe_idx]
                        order = orders[closest_stripe_idx]
                        log.debug(f'fiber {fiber}, order {order} found')

                        fiber = f'fiber_{fiber}'
                        order = f'{order}'

                        if fiber in p_id:
                            p_id[fiber].update({order: p.coefficients})
                        else:
                            p_id[fiber] = {order: p.coefficients}
                    else:
                        log.warning(f'Stripe at {observed_y} could not be identified unambiguously')
                else:
                    log.warning(f'Stripe at {observed_y} could not be identified')

            if debug_level > 1:
                if positions is not None:
                    plt.figure()
                    plt.title("Stripe positions from configuration file")
                    plt.imshow(ad.data[0], origin='lower', vmin=np.min(ad.data[0]), vmax=0.5 * np.max(ad.data[0]))
                    plt.plot(np.repeat(nx / 2, len(y_positions)), y_positions, 'b+', label='original stripe positions')
                    plt.plot(np.repeat(nx / 2, len(y_positions)), y_positions + shift_calculated, 'g+',
                             label=f'corrected stripe positions ({shift_calculated:+.2f})')
                    plt.legend()
                    plt.show()

            if debug_level > 0:
                plt.figure()
                plt.title("Identified orders and fibers")
                plt.imshow(ad.data[0], origin='lower', vmin=np.min(ad.data[0]), vmax=0.5 * np.max(ad.data[0]))
                for f in p_id.keys():
                    for o, p in p_id[f].items():
                        plt.plot(np.poly1d(p)(np.arange(nx)), 'grey')
                        plt.text(0, np.poly1d(p)(0.), f'{f} {o}', fontsize=7,
                                 bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 2})
                plt.tight_layout()
                plt.show()
            ad[0].STRIPES_ID = p_id
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=False)
            return adinputs

    def defineFlatStripes(self,adinputs=None, slit_height=10, extract=False, debug_level=0,**params):
        """
            Extracts and saves flat field profiles.

            For a given slit_height, this function extracts the flat field stripes, calculates a box extracted spectrum and
            normalizes the flat field to generate a 2D profile that is used in the optimal extraction process.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        for ad in adinputs:
            p_id = ad[0].STRIPES_ID
            # # repackage p_id as dict
            # if len(p_id) == 12:  # can do this smarter if access to fiber key identity is kept at end of identifyStripes
            #     p_id = {'fiber_1': dict((colname, p_id[0:6][colname].data) for colname in p_id.colnames),
            #      'fiber_5': dict((colname, p_id[6:][colname].data) for colname in p_id.colnames)}
            # elif len(p_id) == 18:
            #     p_id = {'fiber_2': dict((colname, p_id[0:6][colname].data) for colname in p_id.colnames),
            #             'fiber_3': dict((colname, p_id[6:12][colname].data) for colname in p_id.colnames),
            #             'fiber_4': dict((colname, p_id[12:][colname].data) for colname in p_id.colnames)}
            # elif len(p_id) == 30:
            #     p_id = {'fiber_1': dict((colname, p_id[0:6][colname].data) for colname in p_id.colnames),
            #             'fiber_2': dict((colname, p_id[6:12][colname].data) for colname in p_id.colnames),
            #             'fiber_3': dict((colname, p_id[12:18][colname].data) for colname in p_id.colnames),
            #             'fiber_4': dict((colname, p_id[18:24][colname].data) for colname in p_id.colnames),
            #             'fiber_5': dict((colname, p_id[24:][colname].data) for colname in p_id.colnames)}

            img = ad[0].data
            ny, nx = img.shape
            xx = np.arange(nx)
            index_fiber = np.zeros_like(img, dtype=np.int8)
            index_order = np.zeros_like(img, dtype=np.int8)
            slit_indices_y = np.arange(-slit_height, slit_height).repeat(nx).reshape((2 * slit_height, nx))
            slit_indices_x = np.tile(np.arange(nx), 2 * slit_height).reshape((2 * slit_height, nx))
            for f in p_id.keys():
                for o, p in p_id[f].items():
                    y = np.poly1d(p)(xx)
                    indices = np.rint(slit_indices_y + y).astype(int)
                    valid_indices = np.logical_and(indices < ny, indices > 0)
                    ff = f
                    if isinstance(ff, str):
                        ff = int(''.join(filter(str.isdigit, ff)))
                    index_fiber[indices[valid_indices], slit_indices_x[valid_indices]] = ff
                    oo = o
                    if isinstance(oo, str):
                        oo = int(''.join(filter(str.isdigit, oo)))
                    index_order[indices[valid_indices], slit_indices_x[valid_indices]] = oo

            # image containing only values within stripes, 0 elsewhere
            cleaned_image = np.where(index_order > 0, img, 0)

            if debug_level > 1:
                fig, ax = plt.subplots(1, 3)
                ax[0].imshow(index_fiber, origin='lower')
                ax[1].imshow(index_order, origin='lower')
                ax[2].imshow(cleaned_image, origin='lower')
                plt.show()
            ad[0].INDEX_FIBER = index_fiber.astype(int)
            ad[0].INDEX_ORDER = index_order.astype(int)
            del ad[0].STRIPES_ID # delete interum information
            del ad[0].STRIPES_LOC  # delete interum information
            if extract:
                fiber_tables = []

                for ifib in sorted(p_id.keys(), key=lambda x:x.lower()):
                    # need to sort the keys alphanumerically to ensure proper downstream management of fibers
                    fiber_tables.append(Table.from_pandas(pd.DataFrame(p_id[ifib])))  # loses fiber key identity
                ad[0].STRIPES_ID = vstack(fiber_tables, metadata_conflicts="silent")
            #     # saving full info in DRAGONS, sparse matrix realizations no longer utilized
            #     # ad[0].STRIPES = self._extract_flat_stripes(img, p_id, slit_height, debug_level)
            #
            ad.update_filename(suffix=params['suffix'], strip=True)
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def removeStrayLight(self, adinputs=None,box_size=20, filter_size=20, debug_level=0, **params):
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        adoutputs = []
        for ad in adinputs:
            adint = copy.deepcopy(ad)
            adint2 = copy.deepcopy(ad)
            adout = copy.deepcopy(ad)
            stripes = ad[0].INDEX_FIBER
            orders = ad[0].INDEX_ORDER
            if 'BLUE' in ad.tags:
                log.fullinfo('correcting straylight on blue chip')
                #adint[0].data[bpm == 0] = np.nan
                # Mask the stripes so we only fit the background
                adint[0].data[stripes > 0] = np.nan

                # Add 5 pix to bottom of stripes for orders < 85 to extend masked region in order to account for aberrations
                mask = stripes * 0
                mask[np.where(np.logical_and(orders >= 90, orders <= 124))] = 1
                mask = np.roll(mask, -5, 0)
                adint[0].data[mask == 1] = np.nan

                # on the blue we also miss order 90 -> shift mask of order 91 down 160 pix
                mask = stripes * 0
                mask[orders == 91] = 1
                mask = np.roll(mask, -160, 0)
                mask[2000:, :] = 0
                adint[0].data[mask == 1] = np.nan

                # on the blue we also miss order 125 and 126 -> shift mask of order 124 up 170 pix
                mask = stripes * 0
                mask[orders >= 123] = 1
                mask = np.roll(mask, 170, 0)
                mask[0:2000, :] = 0
                mask[4000:, 3600:] = 1
                adint[0].data[mask == 1] = np.nan
            else:
                log.fullinfo('correcting straylight on red chip')
                #adint[0].data[bpm == 0] = np.nan
                # Mask the stripes so we only fit the background
                adint[0].data[stripes > 0] = np.nan

                # Add 10 pix to top of stripes for orders < 85 to extend masked region in order to account for aberrations
                mask = stripes * 0
                mask[np.where(np.logical_and(orders >= 67, orders <= 85))] = 1
                mask = np.roll(mask, 10, 0)
                adint[0].data[mask == 1] = np.nan

                # Add 5 pix to top of stripes for orders < 85 to extend masked region in order to account for aberrations
                mask = stripes * 0
                mask[np.where(np.logical_and(orders > 85, orders <= 90))] = 1
                mask = np.roll(mask, 5, 0)
                adint[0].data[mask == 1] = np.nan

                # on the red chip we miss fiber order 67 (lower left half of chip)
                # use order 90, shifted down by -200 pix to mask the flux
                mask = stripes * 0
                mask[orders == 67] = 1
                mask = np.roll(mask, -200, 0)
                mask[2000:, :] = 0
                adint[0].data[mask == 1] = np.nan

                # on the red chip we also miss fiber order 94 and 95 (upper right part of chip)
                # use orders 92 and 93, shifted up by 200 pix to mask the flux
                mask = stripes * 0
                mask[orders > 91] = 1
                mask = np.roll(mask, 200, 0)
                mask[:2000, :] = 0
                adint[0].data[mask == 1] = np.nan

            if debug_level > 0:
                plt.figure()
                plt.title('Raw frame')
                plt.imshow(ad[0].data, origin='lower', vmin=0, vmax=100)

                plt.figure()
                plt.title('Raw frame, orders masked')
                plt.imshow(adint[0].data, origin='lower', vmin=0, vmax=100)
                plt.show()

            bkg = Background2D(adint[0].data, (box_size, box_size), filter_size=(filter_size, filter_size), sigma_clip=SigmaClip(sigma=4.),
                               bkg_estimator=MedianBackground(), exclude_percentile=95)

            if debug_level > 2:
                plt.figure()
                plt.title('Raw frame, orders masked')
                plt.imshow(adint[0].data, origin='lower', vmin=0, vmax=200)

                plt.figure()
                plt.title('Background model')
                plt.imshow(bkg.background, origin='lower', vmin=0, vmax=200)

                plt.figure()
                plt.title('Background subtracted data')
                plt.imshow(adint[0].data - bkg.background, origin='lower', vmin=-30, vmax=30)

                plt.figure()
                plt.title('Number of masked pixels in background mesh')
                plt.imshow(bkg.mesh_nmasked, origin='lower', vmin=0, vmax=box_size ** 2)

            # in the blue, we may overshoot in the region where stripes are missing or orders >121
            # a second round of background fitting for negative and near-negative pixels should correct that

            adint2[0].data = ad[0].data - bkg.background
            median_negative = np.median(adint2[0].data[adint2[0].data < 0])
            log.fullinfo(f'Median sub-zero pixel value after 1st iteration: {median_negative}')
            adint2[0].data[adint2[0].data > 2] = np.nan
            adint2[0].data[adint2[0].data< -50] = np.nan
            bkg2 = Background2D(adint2[0].data, (box_size, box_size), filter_size=(filter_size, filter_size), sigma_clip=SigmaClip(sigma=4.),
                                bkg_estimator=MedianBackground(), exclude_percentile=95)

            if debug_level > 2:
                plt.figure()
                plt.title('Background model, correction step')
                plt.imshow(bkg2.background, origin='lower', vmin=-30, vmax=30)

                plt.figure()
                plt.title('Final background subtracted data')
                plt.imshow(ad[0].data - bkg.background - bkg2.background, origin='lower', vmin=-30, vmax=30)

            adout[0].data = ad[0].data - bkg.background - bkg2.background
            log.fullinfo(f'Truncated median sub-zero pixel value: {np.median(adout[0].data[adout[0].data < 0])}')
            adout[0].data[adout[0].data < 0] = 0.01

            adout.update_filename(suffix=params['suffix'], strip=True)
            adoutputs.append(adout)


        gt.mark_history(adoutputs, primname=self.myself(), keyword=timestamp_key)

        return adoutputs

    def separateFlatStreams(self, adinputs=None, **params):
        """
        This primitive splits the flat data into two streams, the 'DFFFD_flats' stream containing DFFFD flats, and main
        containing FDDDF flats. It also warns if non-flats somehow made it into the primitive
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        # Initialize lists of AstroData objects to be added to the streams
        flat_FDDDF_list = []
        flat_DFFFD_list = []
        mislabeled = []
        for ad in adinputs:
            tags = ad.tags
            if "FLAT" in tags and ad.fiber_setup() == ['Flat', 'Dark', 'Dark', 'Dark', 'Flat']:
                flat_FDDDF_list.append(ad)
                log.fullinfo("FDDDF Flat: {}".format(ad.filename))
            elif "FLAT" in tags and ad.fiber_setup() == ['Dark', 'Flat', 'Flat', 'Flat', 'Dark']:
                flat_DFFFD_list.append(ad)
                log.fullinfo("DFFFD Flat: {}".format(ad.filename))
            else:
                mislabeled.append(ad)
                log.warning("Not registered as Flat: {}".format(ad.filename))
        if not flat_FDDDF_list:
            log.warning("No FDDDF Flats in input list")
        if not flat_DFFFD_list:
            log.warning("No DFFFD Flats in input list")

        self.streams["DFFFD_flats"] = flat_DFFFD_list

        return flat_FDDDF_list

    def combineFlatStreams(self, adinputs=None, source=None, **params):
        """
        This primitive recombines the flat data into one master frame, combining the main stream pre-master
        and the 'source' stream pre-master with a simple max comparison at each pix
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        if source not in self.streams.keys():
            log.info("Stream {} does not exist so nothing to transfer".format(source))
            return adinputs

        source_length = len(self.streams[source])
        adinputs_length = len(adinputs)
        if not (adinputs_length == source_length == 1):
            log.warning("Unexpected stream lengths: {} and {}".
                        format(adinputs_length, source_length))
            return adinputs
        adoutputs = []
        adout = copy.deepcopy(adinputs[0])
        adout[0].data = np.max([adinputs[0].data[0], self.streams[source][0].data[0]], axis=0)
        adoutputs.append(adout) # don't need to copy meta-fiber info here, will rerun id'ing on combined image frame
        return adoutputs

    @staticmethod
    def _has_valid_extensions(ad):
        """ Check that the AD has a valid number of extensions. """

        # this needs to be updated at appropriate. 
        return len(ad) in [1]

    def _get_bpm_filename(self, ad):
        """
        Gets bad pixel mask for input MX science frame.
        Returns
        -------
        str/None: Filename of the appropriate bpms
        """
        log = self.log
        arm = ('b' if 'BLUE' in ad.tags else 'r')
        bpm_dir = os.path.join(os.path.dirname(maroonx_maskdb.__file__), 'BPM')
        db_matches = sorted((k, v) for k, v in maroonx_maskdb.bpm_dict.items() if arm in k)

        # If BPM(s) matched, use the one with the latest version number suffix:
        if db_matches:
            bpm = db_matches[-1][1]
        else:
            log.warning('No BPM found for {}'.format(ad.filename))
            return None

        # Prepend standard path if the filename doesn't start with '/'
        return bpm if bpm.startswith(os.path.sep) else os.path.join(bpm_dir, bpm)