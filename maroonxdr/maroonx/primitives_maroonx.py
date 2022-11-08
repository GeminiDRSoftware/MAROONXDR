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
import astrodata
from astrodata.fits import windowedOp
from scipy.ndimage import median_filter
from astropy.table import Table, vstack
import scipy.sparse as sparse
from geminidr.gemini.primitives_gemini import Gemini
from geminidr.core import CCD, NearIR, primitives_preprocess
from astropy.io import fits
from astropy.stats import SigmaClip

from photutils import Background2D, MedianBackground
from .lookups import timestamp_keywords as maroonx_stamps
from .lookups import siddb as maroonx_siddb
from .lookups import maskdb as maroonx_maskdb
from astropy import table
from copy import deepcopy

from scipy.ndimage import measurements  # this set of imports is just for the DQ
from geminidr.gemini.lookups import DQ_definitions as DQ
from astrodata.provenance import add_provenance
from recipe_system.utils.md5 import md5sum

from gempy.gemini import gemini_tools as gt
from gempy.library.nddops import NDStacker
#from gempy.library import astromodels
from . import parameters_maroonx
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

    def addDQ(self, adinputs=None, **params):  # just edited for bpm lookup, can be removed when MX is caldb compliant
        """
        This primitive is used to add a DQ extension to the input AstroData
        object. The value of a pixel in the DQ extension will be the sum of the
        following: (0=good, 1=bad pixel (found in bad pixel mask), 2=pixel is
        in the non-linear regime, 4=pixel is saturated). This primitive will
        trim the BPM to match the input AstroData object(s).

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        static_bpm: str
            Name of bad pixel mask ("default" -> use default from look-up table)
            If set to None, no static_bpm will be added.
        user_bpm: str
            Name of the bad pixel mask created by the user from flats and
            darks.  It is an optional BPM that can be added to the static one.
        illum_mask: bool
            add illumination mask?
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys["addDQ"]
        sfx = params["suffix"]

        # Getting all the filenames first prevents reopening the same file
        # for each science AD
        static_bpm_list = params['static_bpm']
        user_bpm_list = params['user_bpm']

        if static_bpm_list == "default":
            try:
                static_bpm_list = self.caldb.get_processed_bpm(adinputs)
            except:
                static_bpm_list = None
            if static_bpm_list is not None and not all([f is None for f in static_bpm_list.files]):
                static_bpm_list = static_bpm_list.files
            else:
                # TODO once we fully migrate to caldb/server managed bpms, use 2nd line
                # TODO also remove all() check in if above at that time
                static_bpm_list = [self._get_bpm_filename(ad) for ad in adinputs]
                #static_bpm_list = [None] * len(adinputs)

        for ad, static, user in zip(*gt.make_lists(adinputs, static_bpm_list,
                                                   user_bpm_list, force_ad=True)):
            if ad.phu.get(timestamp_key):
                log.warning('No changes will be made to {}, since it has '
                            'already been processed by addDQ'.format(ad.filename))
                continue

            if static is None:
                # So it can be zipped with the AD
                final_static = [None] * len(ad)
            else:
                log.stdinfo("Using {} as static BPM\n".format(static.filename))
                final_static = gt.clip_auxiliary_data(ad, aux=static,
                                                      aux_type='bpm',
                                                      return_dtype=DQ.datatype)

            if user is None:
                final_user = [None] * len(ad)
            else:
                log.stdinfo("Using {} as user BPM".format(user.filename))
                final_user = gt.clip_auxiliary_data(ad, aux=user,
                                                    aux_type='bpm',
                                                    return_dtype=DQ.datatype)

            if static is None and user is None:
                log.stdinfo(f"No BPMs found for {ad.filename} and none supplied by the user.\n")

            for ext, static_ext, user_ext in zip(ad, final_static, final_user):
                if ext.mask is not None:
                    log.warning(f'A mask already exists in extension {ext.id}')
                    continue

                non_linear_level = ext.non_linear_level()
                saturation_level = ext.saturation_level()

                # Need to create the array first for 3D raw F2 data, with 2D BPM
                ext.mask = np.zeros_like(ext.data, dtype=DQ.datatype)
                if static_ext is not None:
                    ext.mask |= static_ext.data
                if user_ext is not None:
                    ext.mask |= user_ext.data

                if saturation_level:
                    log.fullinfo('Flagging saturated pixels in {} extension '
                                 '{} above level {:.2f}'.format(
                                     ad.filename, ext.id, saturation_level))
                    ext.mask |= np.where(ext.data >= saturation_level,
                                         DQ.saturated, 0).astype(DQ.datatype)

                if non_linear_level:
                    if saturation_level:
                        if saturation_level > non_linear_level:
                            log.fullinfo('Flagging non-linear pixels in {} '
                                         'extension {} above level {:.2f}'
                                         .format(ad.filename, ext.id,
                                                 non_linear_level))
                            ext.mask |= np.where((ext.data >= non_linear_level) &
                                                 (ext.data < saturation_level),
                                                 DQ.non_linear, 0).astype(DQ.datatype)
                            # Readout modes of IR detectors can result in
                            # saturated pixels having values below the
                            # saturation level. Flag those. Assume we have an
                            # IR detector here because both non-linear and
                            # saturation levels are defined and nonlin<sat
                            regions, nregions = measurements.label(
                                                ext.data < non_linear_level)
                            # In all my tests, region 1 has been the majority
                            # of the image; however, I cannot guarantee that
                            # this is always the case and therefore we should
                            # check the size of each region
                            region_sizes = measurements.labeled_comprehension(
                                ext.data, regions, np.arange(1, nregions+1),
                                len, int, 0)
                            # First, assume all regions are saturated, and
                            # remove any very large ones. This is much faster
                            # than progressively adding each region to DQ
                            hidden_saturation_array = np.where(regions > 0,
                                                    4, 0).astype(DQ.datatype)
                            for region in range(1, nregions+1):
                                # Limit of 10000 pixels for a hole is a bit arbitrary
                                if region_sizes[region-1] > 10000:
                                    hidden_saturation_array[regions==region] = 0
                            ext.mask |= hidden_saturation_array

                        elif saturation_level < non_linear_level:
                            log.warning(f'{ad.filename} extension {ext.id} '
                                        'has saturation level less than '
                                        'non-linear level')
                        else:
                            log.fullinfo('Saturation and non-linear levels '
                                         'are the same for {}:{}. Only '
                                         'flagging saturated pixels'
                                         .format(ad.filename, ext.id))
                    else:
                        log.fullinfo('Flagging non-linear pixels in {}:{} '
                                     'above level {:.2f}'
                                     .format(ad.filename, ext.id,
                                             non_linear_level))
                        ext.mask |= np.where(ext.data >= non_linear_level,
                                             DQ.non_linear, 0).astype(DQ.datatype)
            if static and static.filename:
                add_provenance(ad, static.filename, md5sum(static.path) or "", self.myself())
            if user and user.filename:
                add_provenance(ad, user.filename, md5sum(user.path) or "", self.myself())

        # Handle latency if reqested
        if params.get("latency", False):
            try:
                adinputs = self.addLatencyToDQ(adinputs, time=params["time"],
                                               non_linear=params["non_linear"])
            except AttributeError:
                log.warning("addLatencyToDQ() not defined in primitivesClass "
                            + self.__class__.__name__)

        # Add the illumination mask if requested
        if params['add_illum_mask']:
            adinputs = self.addIllumMaskToDQ(
                adinputs, **self._inherit_params(params, "addIllumMaskToDQ"))

        # Timestamp and update filenames
        for ad in adinputs:
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=sfx, strip=True)

        return adinputs
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

    def correctImageOrientation(self, adinputs=None, debug_level=0, **params):
        """
        Correct image orientation to proper echelle format.

        flips image so that left lower corner is bluest wavelength, upper right corner is reddest wavelength.
        Echelle orders go from left to right.

        Args:
            adinputs
            debug_level (int): debug level
        Returns:
            adoutputs with correct orientation to nndata
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
                adout[0].mask = np.fliplr(np.flipud(ad[0].mask))

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

    def stackFramesMXCal(self, adinputs=None, **params):
        """
        MX-specific version of stackFrames for calibration frames - changes scaling to average full frame mean to
        purposely scale by etalon flux and its drift between calibration exposures, this function should not be used to
        combine MX science frames.

        Parameters
        ----------
        adinputs : list of :class:`~astrodata.AstroData`
            Any set of 2D.

        suffix : str
            Suffix to be added to output files.

        apply_dq : bool
            Apply DQ mask to data before combining?

        nlow, nhigh : int
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

        operation : str
            Combine method.

        reject_method : str
            Pixel rejection method (none, minmax, sigclip, varclip).

        zero : bool
            Apply zero-level offset to match background levels?

        scale : bool
            Scale images to the same intensity?

        memory : float or None
            Available memory (in GB) for stacking calculations.

        statsec : str
            Section for statistics.

        separate_ext : bool
            Handle extensions separately?

        Returns
        -------
        list of :class:`~astrodata.AstroData`
            Sky stacked image. This list contains only one element. The list
            format is maintained so this primitive is consistent with all the
            others.

        Raises
        ------
        IOError
            If the number of extensions in any of the `AstroData` objects is
            different.

        IOError
            If the shape of any extension in any `AstroData` object is different.

        AssertError
            If any of the `.gain()` descriptors is None.

        AssertError
            If any of the `.read_noise()` descriptors is None.
        """
        def flatten(*args):
            return (el for item in args for el in (
                flatten(*item) if isinstance(item, (list, tuple)) else (item,)))

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys["stackFrames"]
        sfx = params["suffix"]
        memory = params["memory"]
        if memory is not None:
            memory = int(memory * 1000000000)

        zero = params["zero"]
        scale = params["scale"]
        apply_dq = params["apply_dq"]
        separate_ext = params["separate_ext"]
        statsec = params["statsec"]
        reject_method = params["reject_method"]
        save_rejection_map = params["save_rejection_map"]

        if statsec:
            statsec = tuple([slice(int(start)-1, int(end))
                             for x in reversed(statsec.strip('[]').split(','))
                             for start, end in [x.split(':')]])

        if len(adinputs) <= 1:
            log.stdinfo("No stacking will be performed, since at least two "
                        "input AstroData objects are required for stackFrames")
            return adinputs

        if (reject_method == "minmax" and self.mode == "qa" and
                params["nlow"] + params["nhigh"] >= len(adinputs)):
            log.warning("Trying to reject too many images. Setting nlow=nhigh=0.")
            params["nlow"] = 0
            params["nhigh"] = 0

        if len({len(ad) for ad in adinputs}) > 1:
            raise OSError("Not all inputs have the same number of extensions")
        if len({ext.nddata.shape for ad in adinputs for ext in ad}) > 1:
            raise OSError("Not all inputs images have the same shape")

        # We will determine the average gain from the input AstroData
        # objects and add in quadrature the read noise
        gain_list = [ad.gain() for ad in adinputs]
        rn_list = [ad.read_noise() for ad in adinputs]

        # Determine whether we can construct these averages
        process_gain = not None in flatten(gain_list)
        process_rn = not None in flatten(rn_list)

        # Compute gain and read noise of final stacked images
        num_img = len(adinputs)
        num_ext = len(adinputs[0])
        zero_offsets = np.zeros((num_ext, num_img), dtype=np.float32)
        scale_factors = np.ones_like(zero_offsets)

        # Try to determine how much memory we're going to need to stack and
        # whether it's necessary to flush pixel data to disk first
        # Also determine kernel size from offered memory and bytes per pixel
        bytes_per_ext = []
        for ext in adinputs[0]:
            bytes = 0
            # Count _data twice to handle temporary arrays
            bytes += 2 * ext.data.dtype.itemsize
            if ext.variance is not None:
                bytes += ext.variance.dtype.itemsize

            bytes += 2  # mask always created
            bytes_per_ext.append(bytes * np.product(ext.shape))

        if memory is not None and (num_img * max(bytes_per_ext) > memory):
            adinputs = self.flushPixels(adinputs)

        # Compute the scale and offset values by accessing the memmapped data
        # so we can pass those to the stacking function
        # TODO: Should probably be done better to consider only the overlap
        # regions between frames
        if scale or zero:
            levels = np.empty((num_img, num_ext), dtype=np.float32)
            for i, ad in enumerate(adinputs):
                for index in range(num_ext):
                    nddata = (ad[index].nddata.window[:] if statsec is None
                              else ad[index].nddata.window[statsec])
                    scale_mask = (ad[index].nddata.window[:].mask if statsec is None
                              else ad[index].nddata.window[statsec].mask) == DQ.good  # MX specific added line
                    levels[i, index] = np.nansum(nddata.data[scale_mask])  # MX specific changed line
                    # levels[i, index] = gt.measure_bg_from_image(nddata, value_only=True)

            if scale and zero:
                log.warning("Both scale and zero are set. Setting scale=False.")
                scale = False
            if separate_ext:
                # Target value is corresponding extension of first image
                if scale:
                    scale_factors = (np.mean(levels) / levels).T  # MX specific changed line
                    # scale_factors = (levels[0] / levels).T
                else:  # zero=True
                    zero_offsets = (levels[0] - levels).T
            else:
                # Target value is mean of all extensions of first image
                target = np.mean(levels[0])
                if scale:
                    scale_factors = np.tile(target / np.mean(levels, axis=1),
                                            num_ext).reshape(num_ext, num_img)
                else:  # zero=True
                    zero_offsets = np.tile(target - np.mean(levels, axis=1),
                                           num_ext).reshape(num_ext, num_img)
            if scale and np.min(scale_factors) < 0:
                log.warning("Some scale factors are negative. Not scaling.")
                scale_factors = np.ones_like(scale_factors)
                scale = False
            if scale and np.any(np.isinf(scale_factors)):
                log.warning("Some scale factors are infinite. Not scaling.")
                scale_factors = np.ones_like(scale_factors)
                scale = False
            if scale and np.any(np.isnan(scale_factors)):
                log.warning("Some scale factors are undefined. Not scaling.")
                scale_factors = np.ones_like(scale_factors)
                scale = False

        if reject_method == "varclip" and any(ext.variance is None
                                              for ad in adinputs for ext in ad):
            log.warning("Rejection method 'varclip' has been chosen but some"
                        " extensions have no variance. 'sigclip' will be used"
                        " instead.")
            reject_method = "sigclip"

        log.stdinfo("Combining {} inputs with {} and {} rejection"
                    .format(num_img, params["operation"], reject_method))

        stack_function = NDStacker(combine=params["operation"],
                                   reject=reject_method,
                                   log=self.log, **params)

        # NDStacker uses DQ if it exists; if we don't want that, delete the DQs!
        if not apply_dq:
            [setattr(ext, 'mask', None) for ad in adinputs for ext in ad]

        ad_out = astrodata.create(adinputs[0].phu)
        for index, (ext, sfactors, zfactors) in enumerate(
                zip(adinputs[0], scale_factors, zero_offsets)):
            status = (f"Combining extension {ext.id}." if num_ext > 1 else
                      "Combining images.")
            if scale:
                status += " Applying scale factors."
                numbers = sfactors
            elif zero:
                status += " Applying offsets."
                numbers = zfactors
            log.stdinfo(status)
            if (scale or zero) and (index == 0 or separate_ext):
                for ad, value in zip(adinputs, numbers):
                    log.stdinfo(f"{ad.filename:40s}{value:10.3f}")

            shape = adinputs[0][index].nddata.shape
            if memory is None:
                kernel = shape
            else:
                # Chop the image horizontally into equal-sized chunks to process
                # This uses the minimum number of steps and uses minimum memory
                # per step.
                oversubscription = (bytes_per_ext[index] * num_img) // memory + 1
                kernel = ((shape[0] + oversubscription - 1) // oversubscription,) + shape[1:]

            with_uncertainty = True  # Since all stacking methods return variance
            with_mask = apply_dq and not any(ad[index].nddata.window[:].mask is None
                                             for ad in adinputs)
            result = windowedOp(stack_function,
                                [ad[index].nddata for ad in adinputs],
                                scale=sfactors,
                                zero=zfactors,
                                kernel=kernel,
                                dtype=np.float32,
                                with_uncertainty=with_uncertainty,
                                with_mask=with_mask,
                                save_rejection_map=save_rejection_map)
            ad_out.append(result)

            if process_gain:
                gains = [g[index] for g in gain_list]
                # If all inputs have the same gain, the output will also have
                # this gain, and since the header has been copied, we don't
                # need to do anything! (Can't use set() if gains are lists.)
                if not all(g == gains[0] for g in gains):
                    log.warning("Not all inputs have the same gain.")
                    try:
                        output_gain = num_img / np.sum([1/g for g in gains])
                    except TypeError:
                        pass
                    else:
                        ad_out[-1].hdr[ad_out._keyword_for("gain")] = output_gain

            if process_rn:
                # Output gets the rms value of the inputs
                rns = [rn[index] for rn in rn_list]
                output_rn = np.sqrt(np.mean([np.square(np.asarray(rn).mean())
                                             for rn in rns]))
                ad_out[-1].hdr[ad_out._keyword_for("read_noise")] = output_rn

            log.stdinfo("")

        # Propagate REFCAT as the union of all input REFCATs
        refcats = [ad.REFCAT for ad in adinputs if hasattr(ad, 'REFCAT')]
        if refcats:
            try:
                out_refcat = table.unique(table.vstack(refcats, metadata_conflicts='silent'),
                                          keys=('RAJ2000', 'DEJ2000'))
            except KeyError:
                pass
            else:
                out_refcat['Id'] = list(range(1, len(out_refcat)+1))
                ad_out.REFCAT = out_refcat

        # Propagate MDF from first input (no checking that they're all the same)
        if hasattr(adinputs[0], 'MDF'):
            ad_out.MDF = deepcopy(adinputs[0].MDF)

        # Set AIRMASS to be the mean of the input values
        try:
            airmass_kw = ad_out._keyword_for('airmass')
            mean_airmass = np.mean([ad.airmass() for ad in adinputs])
        except Exception:  # generic implementation failure (probably non-Gemini)
            pass
        else:
            ad_out.phu.set(airmass_kw, mean_airmass, "Mean airmass for the exposure")

        # Add suffix to datalabel to distinguish from the reference frame
        if sfx[0] == '_':
            extension = sfx.replace('_', '-', 1).upper()
        else:
            extension = '-' + sfx.upper()
        ad_out.phu.set('DATALAB', "{}{}".format(ad_out.data_label(), extension),
                       self.keyword_comments['DATALAB'])

        # Add other keywords to the PHU about the stacking inputs
        ad_out.orig_filename = ad_out.phu.get('ORIGNAME')
        ad_out.phu.set('NCOMBINE', len(adinputs), self.keyword_comments['NCOMBINE'])
        for i, ad in enumerate(adinputs, start=1):
            ad_out.phu.set('IMCMB{:03d}'.format(i), ad.phu.get('ORIGNAME', ad.filename))

        # Timestamp and update filename and prepare to return single output
        gt.mark_history(ad_out, primname=self.myself(), keyword=timestamp_key)
        ad_out.update_filename(suffix=sfx, strip=True)

        return [ad_out]

    def stackDarks(self, adinputs=None, **params):
        """
        MaroonX-specific version of stack darks allowing scaling for etalon intensity drift that is in MX 'darks'
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        if not all('DARK' in dark.tags for dark in adinputs):
            raise ValueError("Not all inputs have DARK tag")

        if not all(dark.exposure_time() == adinputs[0].exposure_time()
                   for dark in adinputs[1:]):
                raise ValueError("Darks are not of equal exposure time")

        stack_params = self._inherit_params(params, "stackFramesMXCal")
        stack_params.update({'zero': False})  # 'scale': False  # changed line
        adinputs = self.stackFramesMXCal(adinputs, **params)
        return adinputs

    def stackFlats(self, adinputs=None, **params):
        """MaroonX-specific version of stack flats to call correct stackframes for the flats."""
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        stack_params = self._inherit_params(params, "stackFramesMXCal")
        stack_params.update({'zero': False})
        adinputs = self.stackFramesMXCal(adinputs, **stack_params)
        return adinputs

    def findStripes(self,adinputs=None, deg_polynomial=5, median_filter=1, gauss_filter_sigma=3.5, min_peak=0.008,
                    debug_level=0,**params):
        """
        Locates and fits stripes in a flat field or science echelle spectrum.
        Starting in the central column, the algorithm identifies peaks and traces each stripe to the edge of the
        detector by following the brightest pixels along each order. It then fits a polynomial to each stripe. To
        improve algorithm stability, the image is first median filtered and then smoothed with a gaussian. It not only
        eliminates noise, but also ensures that the cross-section profile of the flat becomes peaked in the middle,
        which helps to identify the center of each stripe. Choose gauss_filter accordingly. To avoid false positives,
        only peaks above a certain (relative) intensity threshold are used.

        Parameters
        -------
        adinputs : single MX astrodata object, is either a DFFFD flat, FDDDF flat, or combined FFFFF flat
        deg_polynomial (int): degree of the polynomial fit
        median_filter (int): median filter parameter
        gauss_filter_sigma (float): sigma of the gaussian filter used to smooth the image.
        min_peak (float): minimum relative peak height
        debug_level (int): debug level

        Returns
        -------
        adoutput : single MX astrodata object with STRIPES_LOC extension. This extension temporarily holds the
            fits-unsavable fiber information before it is utilized and then removed.


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
        Requires the findStripes primitive to be run prior during recipe so STRIPES_LOC extension exists

        Parameters
        -------
        adinputs : single MX astrodata object, is either a DFFFD flat, FDDDF flat, or combined FFFFF flat
        positions_dir : lookup array of nominal y positions and fiber/order labels. Shape is Nx3, columns are
        [fibers, orders, y_positions], found in lookups/SID
        debug_level (int): debug level
        selected_fibers (list): fibers illuminated in the flat

        Returns
        ------
        adoutput : single MX astrodata object with STRIPES_ID extension. This extension temporarily holds the
            fits-unsavable fiber information before it is utilized and then removed.

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
        Extracts and saves flat field profiles in a fits-acceptable extension.
        Requires the findStripes and identifyStripes primitives to be run prior during recipe so necessary information
        exists in extensions. Will remove previous STRIPES_ID and STRIPES_LOC extensions and replace with INDEX_FIBER
        and INDEX_ORDER extensions.

        For a given slit_height, this function extracts the flat field stripes, calculates a box extracted spectrum and
        normalizes the flat field to generate a 2D profile that is used in the optimal extraction process.

        Parameters
        -------
        adinputs : single MX astrodata object, is either a DFFFD, FDDDF flat, or combined FFFFF flat
        slit_height (int) : half pixel height of box in each dimension to perform box extraction with
        extract (bool) :  if True, will write STRIPES_ID in fits-acceptable format. Utilized in FFFFF_flat
        debug_level :

        Returns
        -------
        adoutput : single MX astrodata object with INDEX_FIBER, INDEX_ORDER extensions and possibly STRIPES_ID extension
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
            del ad[0].STRIPES_ID  # delete interum information
            del ad[0].STRIPES_LOC  # delete interum information
            if extract:
                fiber_tables = []

                for ifib in sorted(p_id.keys(), key=lambda x:x.lower()):
                    # need to sort the keys alphanumerically to ensure proper downstream management of fibers
                    fiber_tables.append(Table.from_pandas(pd.DataFrame(p_id[ifib])))  # loses fiber key identity
                ad[0].STRIPES_ID = vstack(fiber_tables, metadata_conflicts="silent")
            #     # saving full info in DRAGONS, aka sparse matrix realizations, not utilized not fit compatible
            #     # ad[0].STRIPES = self._extract_flat_stripes(img, p_id, slit_height, debug_level)
            #
            ad.update_filename(suffix=params['suffix'], strip=True)
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def removeStrayLight(self, adinputs=None,box_size=20, filter_size=20, debug_level=0, **params):
        """
        Removes stray light from full frame images for more accurate fiber flux accounting
        Requires the defineStripes primitive to be run prior during recipe so necessary information exists in extensions

        Parameters
        --------
        adinputs : single MX astrodata object, is either a DFFFD or FDDDF flat that has not previously had its
            stray light removed
        box_size (int) : pixel height and width of 'mesh_element' used in background identification sub-routine
        filter_size (int) : pixel height and width of window to perform background identification sub-routine
        debug_level

        Returns
        ------
        adoutput : single MX astrodata object with stray light removed from primary extension
        """
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
                adint[0].data[adint[0].mask != DQ.good] = np.nan  # mask actual bad pixels from bpm
                # Mask the stripes so we only fit the background
                adint[0].data[stripes > 0] = np.nan

                # Add 5 pix to bottom of stripes for orders < 85 to extend masked region in order to account for
                # aberrations
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
                adint[0].data[adint[0].mask != DQ.good] = np.nan  # mask actual bad pixels from bpm
                # Mask the stripes so we only fit the background
                adint[0].data[stripes > 0] = np.nan

                # Add 10 pix to top of stripes for orders < 85 to extend masked region in order to account for
                # aberrations
                mask = stripes * 0
                mask[np.where(np.logical_and(orders >= 67, orders <= 85))] = 1
                mask = np.roll(mask, 10, 0)
                adint[0].data[mask == 1] = np.nan

                # Add 5 pix to top of stripes for orders < 85 to extend masked region in order to account for
                # aberrations
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

            bkg = Background2D(adint[0].data, (box_size, box_size), filter_size=(filter_size, filter_size),
                               sigma_clip=SigmaClip(sigma=4.), bkg_estimator=MedianBackground(), exclude_percentile=95)

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
            adint2[0].data[adint2[0].data < -50] = np.nan
            bkg2 = Background2D(adint2[0].data, (box_size, box_size), filter_size=(filter_size, filter_size),
                                sigma_clip=SigmaClip(sigma=4.), bkg_estimator=MedianBackground(), exclude_percentile=95)

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
        containing FDDDF flats. It also warns if non-flats somehow made it into the list of inputs

        Parameters
        -------
        adinputs : list of MX flats

        Returns
        -------
        'DFFFD_flats' stream
        'main' stream
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

        Parameters
        ------
        'DFFFD_flats' stream : single MX astrodata object
        'main' stream : single MX astrodata object

        Returns
        -------
        adouput : single MX astrodata object with primary extension data from by-pixel-max between stream inputs
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
        adoutputs.append(adout)  # don't need to copy meta-fiber info here, will rerun id'ing on combined image frame
        return adoutputs

    def _get_sid_filename(self, ad):
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

    def _get_bpm_filename(self, ad):
        """
        Gets bad pixel mask for input MX science frame. this function can be removed when MX is bpm caldb compliant
        Returns
        -------
        str/None: Filename of the appropriate bpms
        """
        log = self.log
        arm = ('b' if 'BLUE' in ad.tags else 'r')
        bpm_dir = os.path.join(os.path.dirname(maroonx_maskdb.__file__), 'BPM')
        # db_matches = sorted((k, v) for k, v in maroonx_maskdb.bpm_dict.items() if arm in k)
        #
        # # If BPM(s) matched, use the one with the latest version number suffix:
        # if db_matches:
        #     bpm = db_matches[-1][1]
        # else:
        #     log.warning('No BPM found for {}'.format(ad.filename))
        #     return None
        # Prepend standard path if the filename doesn't start with '/'
        bpm = 'BPM_'+arm+'_0000.fits'
        return bpm if bpm.startswith(os.path.sep) else os.path.join(bpm_dir, bpm)

    @staticmethod
    def _has_valid_extensions(ad):
        """ Check that the AD has a valid number of extensions. """

        # this needs to be updated at appropriate.
        return len(ad) in [1]