#
#                                                                       DRAGONS
#
#                                                         primitives_maroonx.py
# ------------------------------------------------------------------------------
import os

from copy import deepcopy
import pandas as pd
import numpy as np


import astrodata
from astrodata.provenance import add_provenance
from astrodata.fits import windowedOp

from astropy.io import fits
from astropy.stats import SigmaClip
from astropy import table
from astropy.table import Table, vstack

from gempy.gemini import gemini_tools as gt
from gempy.library.nddops import NDStacker
from geminidr.gemini.lookups import DQ_definitions as DQ
from geminidr.gemini.primitives_gemini import Gemini
from geminidr.core import CCD, NearIR

from photutils.background import Background2D, MedianBackground

from scipy.ndimage import median_filter, gaussian_filter
from scipy.ndimage import measurements

from recipe_system.utils.md5 import md5sum
from recipe_system.utils.decorators import parameter_override

import sys
from pathlib import Path
path_root = Path(__file__).parents[1]
sys.path.append(str(path_root))
from maroonx.lookups import timestamp_keywords as maroonx_stamps
from maroonx.lookups import siddb as maroonx_siddb
from maroonx.lookups import maskdb as maroonx_maskdb
from maroonx import parameters_maroonx
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

    def addDQ(self, adinputs=None, **params):
        # just edited for bpm lookup, can be removed when MX is caldb compliant
        """
        This primitive is used to add a DQ extension to the input AstroData
        object. The value of a pixel in the DQ extension will be the sum of the
        following: (0=good, 1=bad pixel (found in bad pixel mask), 2=pixel is
        in the non-linear regime, 4=pixel is saturated). This primitive will
        trim the BPM to match the input AstroData object(s).

        Parameters
        ----------
        adinputs : list of AstroData objects with no DQ extension
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

        Returns
        -------
        adinputs : list of AstroData objects with a DQ extension added to them 
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
                log.warning(f'No changes will be made to {ad.filename}, since it has '
                            'already been processed by addDQ')
                continue

            if static is None:
                # So it can be zipped with the AD
                final_static = [None] * len(ad)
            else:
                log.stdinfo(f"Using {static.filename} as static BPM.\n")
                final_static = gt.clip_auxiliary_data(ad, aux=static,
                                                      aux_type='bpm',
                                                      return_dtype=DQ.datatype)

            if user is None:
                final_user = [None] * len(ad)
            else:
                log.stdinfo(f"Using {user.filename} as user BPM.\n")
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
                    log.fullinfo(f'Flagging saturated pixels in {ad.filename} extension '
                                 f'{ext.id} above level {saturation_level:.2f}')
                    ext.mask |= np.where(ext.data >= saturation_level,
                                         DQ.saturated, 0).astype(DQ.datatype)

                if non_linear_level:
                    if saturation_level:
                        if saturation_level > non_linear_level:
                            log.fullinfo(f'Flagging non-linear pixels in {ad.filename} '
                                         f'extension {ext.id} above level {non_linear_level:.2f}')
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
                                         f'are the same for {ad.filename}:{ext.id}. Only '
                                         'flagging saturated pixels')
                    else:
                        log.fullinfo(f'Flagging non-linear pixels in {ad.filename}:{ext.id} '
                                     f'above level {non_linear_level:.2f}')
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
                log.warning("addLatencyToDQ() not defined in primitivesClass %s"
                            ,self.__class__.__name__)

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
        MAROON-X-specific version of validateData to ignore the invalid WCS
        exception.

        Parameters
        ----------
        adinputs : List of unchecked AstroData objects

        Returns
        -------
        adinputs : List of checked AstroData objects
        """
        try:
            super().validateData(adinputs, suffix=suffix)
        except ValueError as e:
            if 'valid WCS' not in str(e):
                raise
        return adinputs

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

    def checkArm(self, adinputs=None, **params):
        """
        Check that MX frame arm is consistent through all input files, i.e.
        BLUE or RED based on data tags. The first file sets the expected value.
        Currently, assumes 1 astrodata object comes from 1 single-extension
        FITS. Need to update if/when original FITS are MEF.

        Parameters
        ----------
        adinputs - list of un-checked MX frames 

        Returns
        -------
        adoutputs - set of list that passes test,  always at least first frame

        """
        log = self.log
        # find first object's MX-camera
        arm_set = 'BLUE' if 'BLUE' in adinputs[0].tags else \
            'RED' if 'RED' in adinputs[0].tags else 'UNDEFINED'
        if arm_set == 'UNDEFINED':
            log.error(f"{adinputs[0].filename} found without defined camera arm")
            raise IOError()
        adoutputs = []
        if len(adinputs) == 1:
            log.warning('Only one file passed to checkArm')

        # include other objects in list with same tag
        # Warn user and toss frame if not all frames are taken with the same arm
        for ad in adinputs:
            if arm_set not in ad.tags:
                log.warning("Not all frames taken with the same camera arm, "
                            "restricting set to first arm used in list")
                log.warning('tossing frame: %s', ad.filename)
            else:
                ad.update_filename(suffix=params['suffix'], strip=True)
                adoutputs.append(ad)
        return adoutputs

    def correctImageOrientation(self, adinputs=None, **params):
        """
        Correct image orientation to proper echelle format for MAROON-X.
        Flips SCI, if needed, so that left lower corner is bluest wavelength,
        upper right corner is reddest wavelength.
        Resulting echelle orders go from left to right.
        MX blue frames start with incorrect orientation for reduction.
        This primitive must be called after DQ is established and before any
        image arithmetic is performed.

        Parameters
        ----------
        adinputs - list of un-checked MX objects

        Returns
        ----------
        adoutputs - same list as inputs, with correct orientation to SCI
        """

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        adoutputs = []
        for ad in adinputs:
            adout = deepcopy(ad)
            # check for tags inherited from original fits,
            # perform 2-axis flip as needed
            
            #Check if the frame was from the blue arm by looking at the tags
            #If it is, then flip the image
            #TODO: Change this to just look at the arm tag
            if ('BLUE' in ad.tags):
                log.fullinfo(f'{ad.filename} set as blue, orientation flipped')
                adout[0].data = np.fliplr(np.flipud(ad[0].data))
                try:
                    adout[0].mask = np.fliplr(np.flipud(ad[0].mask))
                except:
                    log.warning('DQ not found for %s while orienting image', ad.filename)
            
            #If it is not from the blue arm, then check if it is from the red arm by looking 
            #at the image orientation.  Do not flip the image
            elif ('RED' in ad.tags):
                log.fullinfo(f'{ad.filename} set as red, orientation unchanged')
            
            #In any other case, something has gone wrong- return an error
            else:
                log.error(f"{ad.filename} has no defined orientation")
                raise IOError
            adout.update_filename(suffix=params['suffix'], strip=True)
            adoutputs.append(adout)

        gt.mark_history(adoutputs, primname=self.myself(),
                        keyword=timestamp_key)

        return adoutputs
    
    def addVAR(self, adinputs=None, read_noise=True, poisson_noise = True, **params ):
        """
        Calculates the variance based on the read noise for the chip and the poisson noise
        (the variance in this case is just the number of photons for each pixel).
        The variance is then stored as a FITS extension for each file.

        Parameters
        ----------
        adinputs - list of MX objects without variance extensions
        read_noise - boolean, whether to include read noise in variance calculations
        poisson_noise - boolean, whether to include poisson noise in variance calculations

        Returns
        ----------
        adoutputs - list of MX objects with variance extensions
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        adoutputs = []
        for ad in adinputs:
            log.stdinfo(f'Adding variance extension for {ad.filename}')
            npix_x, npix_y = ad[0].data.shape
            var = np.zeros((npix_x, npix_y))
            if read_noise:
                """
                Check if the frame was from the blue arm by looking at the tags.
                If it is, then use the read noise for the blue arm, 
                otherwise use the read noise for the red arm.
                """
                if ('BLUE' in ad.tags):
                    rn_electron = ad[0] #read noise in electrons for blue arm
                    log.stdinfo(f'{ad.filename} set as Blue, electron read noise is 2.9')
                    
                elif ('RED' in ad.tags):
                    rn_electron = 3.5 #read noise in electrons for red arm
                    log.stdinfo(f'{ad.filename} set as Red, electron read noise is 3.5')
                else:
                    # We should never reach this point, but just in case
                    log.error(f"{ad.filename} has no defined orientation")
                    raise IOError
                read_noise = rn_electron * ad[0].gain()[0] #convert read noise to data units
                var += read_noise**2
            if poisson_noise:
                """
                The variance due to poisson noise is just the number of photons for each pixel.
                Add this to the variance due to the read noise.
                """
                var += ad[0].data
            ad[0].variance = var
            adoutputs.append(ad)
        gt.mark_history(adoutputs, primname=self.myself(),
                        keyword=timestamp_key)
        return adoutputs

    def checkND(self, adinputs=None, **params):
        """
        Check that the ND filter on the sim cal fiber is consistent through
        all MX-input files i.e. illumination is similar intensity as needed for
        good removal. The first file sets the expected value.

        Parameters
        ----------
        adinputs - list of MX-objects

        Returns
        -------
        adoutputs - adinputs that pass test

        """
        log = self.log

        #Get the simcal ND filter setting for the first file
        check_val = adinputs[0].filter_orientation()['ND']
        adoutputs = []

        #In case we have multiple files, check that they all have the same ND filter setting
        if len(adinputs) > 1:
            for ad in adinputs:
                if check_val != ad.filter_orientation()['ND']:
                    log.warning("Not all frames have the same simcal ND filter "
                                "setting, restricting set to first seen")
                else:
                    ad.update_filename(suffix=params['suffix'], strip=True)
                    adoutputs.append(ad)
       
        #If we only have one file with the correct filter setting, return an error
            if len(adoutputs) == 1:
                log.error("Only first frame found, of given, with its"
                          "simcal ND filter setting")
                raise IOError()
       
        #If we only have one file in total, return a warning
        #
        else:
            log.warning('Only one file passed to checkND')
            return adinputs
        return adoutputs

    def stackFramesMXCal(self, adinputs=None, **params):
        """
        MX-specific version of stackFrames for calibration frames - changes
        scaling to average full frame mean to purposely scale by etalon flux
        and its drift between calibration exposures, this function should
        not be used to combine MX science frames.

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
           
        # Check that the input AstroData objects are compatible
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
                              else ad[index].nddata.window[statsec].mask) == DQ.good
                    # MX specific changed line
                    # uses entire good pixel frame, purposely
                    # including etalon flux, to calculate level as sum
                    levels[i, index] = np.nansum(nddata.data[scale_mask])

            if scale and zero:
                log.warning("Both scale and zero are set. Setting scale=False.")
                scale = False
            if separate_ext:
                # Target value is corresponding extension of first image
                if scale:
                    # MX specific changed line
                    # scale each frame by its fractional change from the average
                    # level as opposed to first frame value
                    scale_factors = (np.mean(levels) / levels).T
                else: 
                    zero_offsets = (levels[0] - levels).T
            else:
                # Target value is mean of all extensions of first image
                target = np.mean(levels[0])
                if scale:
                    scale_factors = np.tile(target / np.mean(levels, axis=1),
                                            num_ext).reshape(num_ext, num_img)
                else:  
                    zero_offsets = np.tile(target - np.mean(levels, axis=1),
                                           num_ext).reshape(num_ext, num_img)
            
            # Check for negative, infinite or undefined scale factors
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

        log.stdinfo(f"Combining {num_img} inputs with {params['operation']} and {reject_method} rejection")

        stack_function = NDStacker(combine=params["operation"],
                                   reject=reject_method,
                                   log=self.log, **params)

        # NDStacker uses DQ if it exists; if we don't want that, delete the DQs!
        if not apply_dq:
            for ad in adinputs:
                for ext in ad:
                    setattr(ext, 'mask', None) # delete mask

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
                    # need one digit beyond 10.3f to see differences
                    log.stdinfo(f"{ad.filename:40s}{value:10.4f}")

            shape = adinputs[0][index].nddata.shape
            if memory is None:
                kernel = shape
            else:
                # Chop the image horizontally into equal-sized chunks to process
                # This uses the minimum number of steps and uses minimum memory
                # per step.
                oversubscription = (bytes_per_ext[index] * num_img) // memory + 1
                kernel = ((shape[0] + oversubscription - 1) // oversubscription,) + shape[1:]

            with_mask = apply_dq and not any(ad[index].nddata.window[:].mask is None
                                             for ad in adinputs)
            result = windowedOp(stack_function,
                                [ad[index].nddata for ad in adinputs],
                                scale=sfactors,
                                zero=zfactors,
                                kernel=kernel,
                                dtype=np.float32,
                                with_uncertainty=True,
                                with_mask=with_mask,
                                save_rejection_map=save_rejection_map)
            ad_out.append(result)
            log.stdinfo("")

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
        ad_out.phu.set('DATALAB', f"{ad_out.data_label()}{extension}",
                       self.keyword_comments['DATALAB'])

        # Add other keywords to the Fits Header Unit about the stacking inputs
        ad_out.orig_filename = ad_out.phu.get('ORIGNAME')
        ad_out.phu.set('NCOMBINE %d %s', len(adinputs), self.keyword_comments['NCOMBINE'])
        for i, ad in enumerate(adinputs, start=1):
            ad_out.phu.set('IMCMB{:03d}'.format(i), ad.phu.get('ORIGNAME', ad.filename))

        # Timestamp and update filename and prepare to return single output
        gt.mark_history(ad_out, primname=self.myself(), keyword=timestamp_key)
        ad_out.update_filename(suffix=sfx, strip=True)

        return [ad_out]

    def stackDarks(self, adinputs=None, **params):
        """
        MX-specific version of stack darks allowing scaling for etalon
        intensity drift that is in MX 'darks'
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        # Check that all inputs are DARKs and have the same exposure time - not MX specific
        if not all('DARK' in dark.tags for dark in adinputs):
            raise ValueError("Not all inputs have DARK tag")

        if not all(dark.exposure_time() == adinputs[0].exposure_time()
                   for dark in adinputs[1:]):
            raise ValueError("Darks are not of equal exposure time")

        # MX specific-changed lines start
        # MX 'dark' frames have flux in them, need to scale.
        # Also utilizes special stackFramesMXCal scaling.
        stack_params = self._inherit_params(params, "stackFramesMXCal")
        stack_params.update({'zero': False})
        adinputs = self.stackFramesMXCal(adinputs, **params)
        # MX specific-changed lines end
        
        return adinputs

    def stackFlats(self, adinputs=None, **params):
        """
        MaroonX-specific version of stack flats to call correct stackframes
        for the flats.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        # MX specific-changed line
        # Utilizes special stackFramesMXCal scaling.
        stack_params = self._inherit_params(params, "stackFramesMXCal")
        stack_params.update({'zero': False})
        adinputs = self.stackFramesMXCal(adinputs, **stack_params)
        return adinputs

    def findStripes(self, adinputs=None, deg_polynomial=5, med_filter=1,
                    gauss_filter_sigma=3.5, min_peak=0.008, **params):
        """
        Locates and fits stripes in a flat field spectrum.
        Starting in the central column, the algorithm identifies peaks and
        traces each stripe to the edge of the detector by following the
        brightest pixels along each order. It then fits a polynomial to each
        stripe. To improve algorithm stability, the image is first median
        filtered and then smoothed with a gaussian. It not only eliminates
        noise, but also ensures that the cross-section profile of the flat
        becomes peaked in the middle, which helps to identify the center of each
        stripe. Choose gauss_filter accordingly. To avoid false positives, only
        peaks above a certain (relative) intensity threshold are used.

        Parameters
        -------
        adinputs : single MX astrodata object, is either a DFFFD flat,
            FDDDF flat, or combined FFFFF flat
        deg_polynomial : degree of the polynomial fit
        med_filter : median filter parameter
        gauss_filter_sigma : sigma of the gaussian filter used to
            smooth the image.
        min_peak : minimum relative peak height

        Returns
        -------
        adoutput : single MX astrodata object with STRIPES_LOC extension.
        This extension temporarily holds the fits-unsavable fiber information
        before it is utilized and then removed.


        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:

            npix_y, npix_x = ad.data[0].shape
            # first smooth frame to improve algorithm stability
            filtered_ad = median_filter(ad.data[0], med_filter)
            filtered_ad = gaussian_filter(filtered_ad, gauss_filter_sigma)

            # find peaks in center column
            data = filtered_ad[:, int(npix_x / 2)]
            peaks = np.r_[True, data[1:] >= data[:-1]] & \
                    np.r_[data[:-1] > data[1:], True]

            idx = np.logical_and(peaks, data > min_peak * np.max(data))
            maxima = np.arange(npix_y)[idx]

            # filter out maxima too close to the boundary to avoid problems
            maxima = maxima[maxima > 3]
            maxima = maxima[maxima < npix_y - 3]

            n_order = len(maxima)
            log.fullinfo(f'Number of stripes found: {n_order}')

            orders = np.zeros((n_order, npix_x))
            # walk through to the left and right along the maximum of the order
            for row_max, row in enumerate(maxima):
                column = int(npix_x / 2)
                orders[row_max, column] = row
                start_row = row
                # walk right
                while column + 1 < npix_x:
                    column += 1
                    args = [start_row]
                    if start_row - 1 > 1:
                        args.append(start_row - 1)
                    else:
                        args.append(1)
                    if start_row + 1 < npix_y - 1:
                        args.append(start_row + 1)
                    else:
                        args.append(npix_y - 1)
                    pixel = filtered_ad[args, column]
                    # new maximum
                    start_row = args[int(np.argmax(pixel))]
                    orders[row_max, column] = start_row

                column = int(npix_x / 2)
                start_row = row
                # walk left
                while column > 0:
                    column -= 1
                    args = [start_row]
                    if start_row - 1 > 1:
                        args.append(start_row - 1)
                    else:
                        args.append(1)
                    if start_row + 1 < npix_y - 1:
                        args.append(start_row + 1)
                    else:
                        args.append(npix_y - 1)
                    pixel = filtered_ad[args, column]
                    # new maximum
                    start_row = args[int(np.argmax(pixel))]
                    orders[row_max, column] = start_row

            # do Polynomial fit for each order
            log.fullinfo(f'Fit polynomial of order {deg_polynomial} '
                         f'to each stripe')
            x_pixels = np.arange(npix_x)
            polynomials = [np.poly1d(np.polyfit(x_pixels, o, deg_polynomial))
                           for o in orders]
            # polynomials is FITs-unsaveable dict of dicts,
            # will use and then remove before storing
            ad[0].STRIPES_LOC = polynomials
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=False)
        return adinputs

    def identifyStripes(self, adinputs=None, positions_dir=None,
                        selected_fibers=None, **params):
        """
        Identifies the stripes by assigning their proper order and fiber number,
        including correction for the possibilitiy that the spectra have shifted
        up/down in the cross-dispersion direction since the reference was made.
        Requires the findStripes primitive to be run prior during recipe so
        the stripes are located in the input, i.e. STRIPES_LOC extension exists.

        Parameters
        -------
        adinputs : single MX astrodata object, is either a DFFFD flat,
            FDDDF flat, or combined FFFFF flat with STRIPES_LOC extension
        positions_dir : lookup fits location of nominal y positions and
            fiber/order labels. Shape is Nx3, columns are
            [fibers, orders, y_positions], nominally found in lookups/SID
        selected_fibers : fibers illuminated in the flat, if None assumes all
            can work if not given on partially illuminated frame, but best
            practice is to explicitly identify on function call.

        Returns
        ------
        adoutput : single MX astrodata object with STRIPES_ID extension.
        This extension temporarily holds the fits-unsavable fiber information
        before it is utilized and then removed. A new extension REMOVED_STRIPES
        also saves the polynomial info for every stripe that is not identified
        from the original set inherited from findStripes

        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            # first obtain reference positions (position_dir)
            if positions_dir is None:
                positions_dir = self._get_sid_filename(ad)
                log.info(positions_dir)
            with fits.open(positions_dir, 'readonly') as lookup_frame:
                positions = lookup_frame[1].data
            _, npix_x = ad.data[0].shape
            p_id = {}

            # Create the list of illuminated fibers from a comma separated string
            selected_fibers = list(np.asarray((selected_fibers.split(','))
                                              ).astype(int))
            
            # Boolean to check if all fibers are being used (selected_fibers=None)
            use_all_fibers = selected_fibers is None

            y_positions = positions['fiber_position'].astype('int')
            fibers = positions['identify_fiber'].astype('int')
            orders = positions['fiber_order'].astype('int')
            unique_fibers = np.unique(fibers)

            # If the selected fibers are not given, use all fibers
            if selected_fibers is None:
                selected_fibers = unique_fibers

            # Otherwise, remove the none illuminated fibers
            for fiber_iter in unique_fibers:
                if fiber_iter not in selected_fibers:
                    idx = np.where(fibers == fiber_iter)
                    y_positions = np.delete(y_positions, idx)
                    fibers = np.delete(fibers, idx)
                    orders = np.delete(orders, idx)

            used = np.zeros_like(y_positions)
            # second, pull in previously calculated locations on frame
            observed_y = []
            for _, poly in enumerate(ad[0].STRIPES_LOC):
                observed_y.append(np.poly1d(poly)(npix_x / 2))

            # third, compare the reference values plus an arbitrary shift
            # to the on frame location, over a range in shifts
            shifts = np.linspace(-200, 200, 2000)
            total_distances = []
            for shift in shifts:
                distances = []
                for y in observed_y:
                    closest_stripe_idx = np.argmin(np.abs(y_positions +
                                                          shift - y))
                    distances.append(np.abs(y_positions[closest_stripe_idx] +
                                            shift - y))
                total_distances.append(np.array(distances).sum())

            # very important correction for fiber IDs on subsets of the
            # illuminated fibers, weights towards minima closer to the
            # initial guess
            if not use_all_fibers:
                total_distances += np.abs(shifts) * 2
            # fourth, pick the best-aligned arbitrary shift reference for the
            # all stripes
            shift_calculated = shifts[np.argmin(total_distances)]
            # finally, create dict of dicts assigning stripe fibers and orders
            # to those found in the flat frame.
            # Additionally, save unidentified stripe polynomial info in
            # FITs savable REMOVED_STRIPES extension
            unidentified_stripes = []
            for _, poly in enumerate(ad[0].STRIPES_LOC): # for each stripe
                # Observed y position of stripe is given by the stored polynomial and the  center column of frame
                observed_y = np.poly1d(poly)(npix_x / 2)
                # Find the closest stripe in the reference frame, using the SID file's y positions
                closest_stripe_idx = np.argmin(np.abs(y_positions +
                                                      shift_calculated -
                                                      observed_y))
                
                # If the stripe is within 7 pixels of the reference, assign it
                if np.abs(y_positions[closest_stripe_idx] + shift_calculated -
                          observed_y) < 7:
                    if used[closest_stripe_idx] == 0:
                        used[closest_stripe_idx] = 1
                        fiber = fibers[closest_stripe_idx]
                        order = orders[closest_stripe_idx]
                        log.debug(f'fiber {fiber}, order {order} found')

                        fiber = f'fiber_{fiber}'
                        order = f'{order}'

                        if fiber in p_id:
                            p_id[fiber].update({order: poly.coefficients})
                        else:
                            p_id[fiber] = {order: poly.coefficients}
                    # If the stripe is not in the dictionary, warn the user
                    else:
                        log.warning(f'Stripe at {observed_y} could not be identified unambiguously')
                        unidentified_stripes.append(poly.coefficients)
                # If the stripe is not within 7 pixels of the reference, warn the user, and add it to the unidentified stripes
                else:
                    log.warning(f'Stripe at {observed_y} could not be identified')
                    unidentified_stripes.append(poly.coefficients)
            # p_id is FITs-unsavable dict of dicts,
            # will use and then reformat before storing
            ad[0].STRIPES_ID = p_id
            # REMOVED_STRIPES is FITS-savable np.array of unidentified stripes
            # there is no meta-info being dropped as the stripe is unidentified
            ad[0].REMOVED_STRIPES = np.array(unidentified_stripes)
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=False)
            return adinputs

    def defineFlatStripes(self,adinputs=None, slit_height=10, extract=False,
                          **params):
        """
        Saves fiber location info based on flat field info for stray light
        removal (extract=False) and for future science extraction (extract=True).
        Requires the findStripes and identifyStripes primitives to be run prior
        during recipe so necessary information exists in input extensions.
        Will remove previous (improperly formatted, but fast)
        STRIPES_ID and STRIPES_LOC extensions and replace with INDEX_FIBER
        and INDEX_ORDER pixel map extensions, as needed in straylight removal,
        and (if extract=True) a FITS savable STRIPES_ID and STRIPES_FIBERS.

        For a given slit_height, this function extracts the flat field stripes,
        calculates a box extracted spectrum and normalizes the flat field to
        generate a 2D pixel map that is used in the straylight removal.

        STRIPES_ID and STRIPES_FIBERS contain the by-spectral-order
        polynomial plate solution for each illuminated fiber that is utilized
        to define 2D extraction regions in science extractions.

        Parameters
        -------
        adinputs : single MX astrodata object, is either a DFFFD, FDDDF flat,
            or combined FFFFF flat
        slit_height : half pixel height of box in each dimension to
            perform box extraction with
        extract :  if True, will write STRIPES_ID in fits-acceptable
            format. Utilized in combined, all fiber illuminated FFFFF_flat

        Returns
        -------
        adoutput : single MX astrodata object with INDEX_FIBER, INDEX_ORDER
            extensions and possibly STRIPES_ID and STRIPES_FIBERS extensions
        """

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        for ad in adinputs:
            p_id = ad[0].STRIPES_ID #Create dictionary of stripes
            img = ad[0].data 
            npix_y, npix_x = img.shape
            x_pixels = np.arange(npix_x)
            index_fiber = np.zeros_like(img, dtype=np.int8)  # 2D pixel map
            index_order = np.zeros_like(img, dtype=np.int8) 
            slit_indices_y = np.arange(-slit_height, slit_height)\
                .repeat(npix_x).reshape((2 * slit_height, npix_x)) 
            slit_indices_x = np.tile(np.arange(npix_x), 2 * slit_height)\
                .reshape((2 * slit_height, npix_x))
            for fiber_iter in p_id.keys():
                for order_iter, poly in p_id[fiber_iter].items():
                    y = np.poly1d(poly)(x_pixels)
                    indices = np.rint(slit_indices_y + y).astype(int)
                    valid_indices = np.logical_and(indices < npix_y,
                                                   indices > 0)
                    final_fiber = fiber_iter
                    if isinstance(final_fiber, str):
                        final_fiber = int(''.join(filter(str.isdigit,
                                                         final_fiber)))
                    index_fiber[indices[valid_indices],
                                slit_indices_x[valid_indices]] = final_fiber
                    final_order = order_iter
                    if isinstance(final_order, str):
                        final_order = int(''.join(filter(str.isdigit,
                                                         final_order)))
                    index_order[indices[valid_indices],
                                slit_indices_x[valid_indices]] = final_order

            ad[0].INDEX_FIBER = index_fiber.astype(int)
            ad[0].INDEX_ORDER = index_order.astype(int)
            del ad[0].STRIPES_ID  # delete interim information in improper format
            del ad[0].STRIPES_LOC  # delete interim information

            if extract:
                fiber_tables = []
                for ifib in sorted(p_id.keys(), key=lambda x: x.lower()):
                    fiber_tables.append(
                        Table.from_pandas(pd.DataFrame(p_id[ifib])))
                ad[0].STRIPES_ID = vstack(fiber_tables,
                                          metadata_conflicts="silent")
                ad[0].STRIPES_FIBERS = np.array([key[-1] for key in list(sorted(p_id.keys(), key=lambda x: x.lower()))]).astype(int)
                # actually extracting stripes using full info in DRAGONS,
                # aka sparse matrix realizations,
                # not utilized as csc not FITS compatible

            ad.update_filename(suffix=params['suffix'], strip=True)
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def removeStrayLight(self, adinputs=None, box_size=21, filter_size=21,
                         snapshot=False, **params):
        """
        Removes stray light from full frame images for more accurate fiber flux
        accounting. Requires the defineStripes primitive to be run prior during
        recipe so INDEX_FIBER and INDEX_ORDER extensions exist to define pixel
        locations across frame within fiber traces to avoid when finding stray
        light.

        Parameters
        --------
        adinputs : single MX astrodata object, is either a DFFFD or FDDDF flat
            that has not previously had its stray light removed
        box_size : pixel height and width of 'mesh_element' used in
            background identification sub-routine
        filter_size : pixel height and width of window to perform
            background identification sub-routine
        snapshot : Bool to save difference frame of removed stray light as
            extension STRAYLIGHT_DIFFERENCE

        Returns
        ------
        adoutput : single MX astrodata object with stray light removed from
            SCI
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        adoutputs = []
        for ad in adinputs:
            adint = deepcopy(ad)
            adint2 = deepcopy(ad)
            adout = deepcopy(ad)
            stripes = ad[0].INDEX_FIBER
            orders = ad[0].INDEX_ORDER

            log.fullinfo('correcting straylight ')
            # Mask all actual bad pixels from bpm
            adint[0].data[adint[0].mask != DQ.good] = np.nan
            # Mask all the stripes so we only fit the background
            adint[0].data[stripes > 0] = np.nan

            #Check if we are on the blue chip
            if 'BLUE' in ad.tags:
                log.fullinfo('on blue chip')
                # Add 5 pix to bottom of stripes for orders < 85 to extend
                # masked region in order to account for aberrations
                mask = stripes * 0
                mask[np.where(np.logical_and(orders >= 90, orders <= 124))] = 1
                mask = np.roll(mask, -5, 0)
                adint[0].data[mask == 1] = np.nan

                # on the blue we also miss order 90 -> shift mask of
                # order 91 down 160 pix
                mask = stripes * 0
                mask[orders == 91] = 1
                mask = np.roll(mask, -160, 0)
                mask[2000:, :] = 0
                adint[0].data[mask == 1] = np.nan

                # on the blue we also miss order 125 and 126 -> shift mask
                # of order 124 up 170 pix
                mask = stripes * 0
                mask[orders >= 123] = 1
                mask = np.roll(mask, 170, 0)
                mask[0:2000, :] = 0
                mask[4000:, 3600:] = 1
                adint[0].data[mask == 1] = np.nan

            #Else we must be on the red chip
            else:
                log.fullinfo('on red chip')
                # Add 10 pix to top of stripes for orders < 85 to extend
                # masked region in order to account for aberrations
                mask = stripes * 0
                mask[np.where(np.logical_and(orders >= 67, orders <= 85))] = 1
                mask = np.roll(mask, 10, 0)
                adint[0].data[mask == 1] = np.nan

                # Add 5 pix to top of stripes for orders < 85 to extend
                # masked region in order to account for aberrations
                mask = stripes * 0
                mask[np.where(np.logical_and(orders > 85, orders <= 90))] = 1
                mask = np.roll(mask, 5, 0)
                adint[0].data[mask == 1] = np.nan

                # on the red chip we miss fiber order 67 (lower left half of
                # chip) use order 90, shifted down by -200 pix to mask the flux
                mask = stripes * 0
                mask[orders == 67] = 1
                mask = np.roll(mask, -200, 0)
                mask[2000:, :] = 0
                adint[0].data[mask == 1] = np.nan

                # on the red chip we also miss fiber order 94 and 95 (upper
                # right part of chip) use orders 92 and 93, shifted up by
                # 200 pix to mask the flux
                mask = stripes * 0
                mask[orders > 91] = 1
                mask = np.roll(mask, 200, 0)
                mask[:2000, :] = 0
                adint[0].data[mask == 1] = np.nan

            # Create a background mesh with box size and filter size 21 
            # (changed from 20 due to odd number requirements) 
            bkg = Background2D(adint[0].data, (box_size, box_size),
                               filter_size=(filter_size, filter_size),
                               sigma_clip=SigmaClip(sigma=4.),
                               bkg_estimator=MedianBackground(),
                               exclude_percentile=95) 

            # in the blue, we may overshoot in the region where stripes are
            # missing or orders >121 a second round of background fitting
            # for negative and near-negative pixels should correct that

            adint2[0].data = ad[0].data - bkg.background
            median_negative = np.median(adint2[0].data[adint2[0].data < 0])
            log.fullinfo(f'Median sub-zero pixel value after 1st iteration: '
                         f'{median_negative}')
            adint2[0].data[adint2[0].data > 2] = np.nan
            adint2[0].data[adint2[0].data < -50] = np.nan
            bkg2 = Background2D(adint2[0].data, (box_size, box_size),
                                filter_size=(filter_size, filter_size),
                                sigma_clip=SigmaClip(sigma=4.),
                                bkg_estimator=MedianBackground(),
                                exclude_percentile=95)
            
            # Perform the second background fitting
            adout[0].data = ad[0].data - bkg.background - bkg2.background 
            log.fullinfo(f'Truncated median sub-zero pixel value: '
                         f'{np.median(adout[0].data[adout[0].data < 0])}')
            adout[0].data[adout[0].data < 0] = 0.01
            if snapshot:
                # used for unit testing, save straylight difference
                # but return original data so it can be compared in test
                adout[0].STRAYLIGHT_DIFFERENCE = adout[0].data - \
                                                 deepcopy(ad)[0].data
                adout[0].data = deepcopy(ad)[0].data
            adout.update_filename(suffix=params['suffix'], strip=True)
            adoutputs.append(adout)

        gt.mark_history(adoutputs, primname=self.myself(),
                        keyword=timestamp_key)

        return adoutputs

    def separateFlatStreams(self, adinputs=None, **params):
        """
        This primitive splits the flat data into two streams, the 'DFFFD_flats'
        stream containing DFFFD flats, and main containing FDDDF flats.
        It also warns if non-flats somehow made it into the list of inputs

        Parameters
        -------
        adinputs : list of MX flats
        **params needed for access to stream
        Returns
        -------
        'DFFFD_flats' stream
        'main' stream
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        # Initialize lists of AstroData objects to be added to the streams
        flat_fdddf_list = []
        flat_dfffd_list = []
        mislabeled = []
        for ad in adinputs:
            tags = ad.tags
            #Create list of FDDDF flats to go in the main stream
            if "FLAT" in tags and ad.fiber_setup() == ['Flat', 'Dark', 'Dark',
                                                       'Dark', 'Flat']:
                flat_fdddf_list.append(ad)
                log.fullinfo(f"FDDDF Flat: {ad.filename}")
            #Create list of DFFFD flats to go in the DFFFD_flats stream
            elif "FLAT" in tags and ad.fiber_setup() == ['Dark', 'Flat', 'Flat',
                                                         'Flat', 'Dark']:
                flat_dfffd_list.append(ad)
                log.fullinfo(f"DFFFD Flat: {ad.filename}")
            #Warn if non-flats are in the input list- any other fiber setup is incorrect
            else:
                mislabeled.append(ad)
                log.warning(f"Not registered as Flat: {ad.filename}")
            #Provide warnings if we do not have both types of flats
        if not flat_fdddf_list:
            log.warning("No FDDDF Flats in input list")
        if not flat_dfffd_list:
            log.warning("No DFFFD Flats in input list")

        self.streams["DFFFD_flats"] = flat_dfffd_list
        return flat_fdddf_list

    def combineFlatStreams(self, adinputs=None, source=None, **params):
        """
        Recombines the flat data into one processed frame,
        combining the main stream pre-processed and the 'source' stream
        pre-processed with a simple max comparison at each pix

        Parameters
        ------
        'DFFFD_flats' stream : single MX astrodata object
        'main' stream : single MX astrodata object
        **params needed for access to stream
        Returns
        -------
        adoutput : single FFFFF_flat MX astrodata object with primary extension
            data as combined all fiber illuminated flat
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting")) #Log the start of the primitive

        if source not in self.streams.keys(): #Check the streams dictionary to see if there is source exists
            log.info("Stream %s does not exist so nothing to "
                     "transfer", source)
            return adinputs #If source does not exist, return the provided input without modification

        source_length = len(self.streams[source]) #Get the length of the source stream
        adinputs_length = len(adinputs) #Get the length of the input stream

        #We expect the source stream to have a length of 1, and the input stream to have a length of 1
        #as we have a single DFFFD flat and a single FDDDF flat. If this is not the case, we log a warning.

        if not adinputs_length == source_length == 1: 
            log.warning("Unexpected stream lengths: %s and %s",
                        adinputs_length, source_length)
            #Return the input without modification as we have unexpected stream lengths
            return adinputs 
        #Provided the stream lengths are as expected, we can proceed with the combination
        adoutputs = []
        adout = deepcopy(adinputs[0]) 
        #Combine the data from the two streams by taking the max at each pixel
        adout[0].data = np.max([adinputs[0].data[0],
                                self.streams[source][0].data[0]], axis=0) 
        # Do the same for the variance
        adout[0].variance = np.max([adinputs[0].variance[0],
                                    self.streams[source][0].variance[0]], axis=0)
        adoutputs.append(adout)
        # For the rest of the extensions, we do not need to do this because we
        # will rerun id'ing on combined image frame
        return adoutputs

    def _get_sid_filename(self, ad):
        """
        Gets stripe ID file for input frame.  SID will not be caldb compliant as it is 
        instrument specific.

        Returns
        -------
        str/None: Filename of appropriate sid
        """
        log = self.log
        arm = ('b' if 'BLUE' in ad.tags else 'r') #Get appropriate arm
        sid_dir = os.path.join(os.path.dirname(maroonx_siddb.__file__), 'SID')
        db_matches = sorted((k, v) for k, v in maroonx_siddb.sid_dict.items()
                            if arm in k) #Check if there is a Stripe ID file for the given arm
        if db_matches:
            sid = db_matches[-1][1]
        else:
            log.warning('No SID found for %s',ad.filename)
            return None
        return sid if sid.startswith(os.path.sep) else \
            os.path.join(sid_dir, sid)

    def _get_bpm_filename(self, ad):
        """
        Gets bad pixel mask for input MX science frame.
        this function can be removed when MX is bpm caldb compliant

        Returns
        -------
        str/None: Filename of the appropriate bpms
        """
        arm = ('b' if 'BLUE' in ad.tags else 'r') #Get appropriate arm
        bpm_dir = os.path.join(os.path.dirname(maroonx_maskdb.__file__), 'BPM')
        bpm = 'BPM_'+arm+'_0000.fits' #Append appropriate arm to the bpm name
        return bpm if bpm.startswith(os.path.sep) else \
            os.path.join(bpm_dir, bpm)

    @staticmethod
    def _has_valid_extensions(ad):
        """ Check that the AD has a valid number of extensions. """

        # this needs to be updated as appropriate.
        return len(ad) in [1]
