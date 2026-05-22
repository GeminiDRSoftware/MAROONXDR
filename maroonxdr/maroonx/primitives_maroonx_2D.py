"""Maroon-X primitives to reduce 2D spectra."""

# ------------------------------------------------------------------------------
from copy import deepcopy

import astrodata
import matplotlib.pyplot as plt
from maroonxdr.maroonx.maroonx_plots import plot_backgroundfit
import numpy as np
import pandas as pd
from astrodata.fits import windowedOp
from astrodata.provenance import add_provenance
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.table import Table, vstack
from geminidr.core import CCD, NearIR
from geminidr.gemini.lookups import DQ_definitions as DQ
from geminidr.gemini.primitives_gemini import Gemini
from gempy.gemini import gemini_tools as gt
from gempy.library.nddops import NDStacker
from matplotlib.backends.backend_pdf import PdfPages
from photutils.background import Background2D, MedianBackground
from recipe_system.utils.decorators import parameter_override
from recipe_system.utils.md5 import md5sum
from scipy.ndimage import gaussian_filter, measurements, median_filter

from . import maroonx_utils, parameters_maroonx_2D
from .lookups import timestamp_keywords as maroonx_stamps
from .primitives_calibdb_maroonx import CalibDBMAROONX


@parameter_override
class MAROONX(CalibDBMAROONX, Gemini, CCD, NearIR):
    """Any primitives specific to MAROON-X can go here."""

    tagset = {'GEMINI', 'MAROONX'}

    def _initialize(self, adinputs, **kwargs):
        self.inst_lookups = 'maroonxdr.maroonx.lookups'
        super()._initialize(adinputs, **kwargs)
        self._param_update(parameters_maroonx_2D)
        self.timestamp_keys.update(maroonx_stamps.timestamp_keys)

    def addDQ(self, adinputs=None, **params):
        # just edited for bpm lookup, can be removed when MX is caldb compliant
        """Add a DQ extension to the input AstroData objects.

        The value of a pixel in the DQ extension will be the sum of the
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
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        timestamp_key = self.timestamp_keys['addDQ']
        sfx = params['suffix']

        # Getting all the filenames first prevents reopening the same file
        # for each science AD
        static_bpm_list = params['static_bpm']
        user_bpm_list = params['user_bpm']

        if static_bpm_list == 'default':
            try:
                static_bpm_list = self.caldb.get_processed_bpm(adinputs)
            except:
                static_bpm_list = None
            if static_bpm_list is not None and not all(
                [f is None for f in static_bpm_list.files]
            ):
                static_bpm_list = static_bpm_list.files
            else:
                # TODO once we fully migrate to caldb/server managed bpms, use 2nd line
                # TODO also remove all() check in if above at that time
                static_bpm_list = [
                    maroonx_utils.get_bpm_filename(ad) for ad in adinputs
                ]
                # static_bpm_list = [None] * len(adinputs) # noqa

        for ad, static, user in zip(
            *gt.make_lists(adinputs, static_bpm_list, user_bpm_list, force_ad=True)
        ):
            if ad.phu.get(timestamp_key):
                log.warning(
                    'No changes will be made to %s, since it has '
                    'already been processed by addDQ',
                    ad.filename,
                )
                continue

            if static is None:
                # So it can be zipped with the AD
                final_static = [None] * len(ad)
            else:
                log.stdinfo(f'Using {static.filename} as static BPM.\n')
                final_static = gt.clip_auxiliary_data(
                    ad, aux=static, aux_type='bpm', return_dtype=DQ.datatype
                )

            if user is None:
                final_user = [None] * len(ad)
            else:
                log.stdinfo(f'Using {user.filename} as user BPM.\n')
                final_user = gt.clip_auxiliary_data(
                    ad, aux=user, aux_type='bpm', return_dtype=DQ.datatype
                )

            if static is None and user is None:
                log.stdinfo(
                    f'No BPMs found for {ad.filename} and none supplied by the user.\n'
                )

            for ext, static_ext, user_ext in zip(ad, final_static, final_user):
                if ext.mask is not None:
                    log.warning('A mask already exists in extension %s', ext.id)
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
                    log.fullinfo(
                        f'Flagging saturated pixels in {ad.filename} extension '
                        f'{ext.id} above level {saturation_level:.2f}'
                    )
                    ext.mask |= np.where(
                        ext.data >= saturation_level, DQ.saturated, 0
                    ).astype(DQ.datatype)

                if non_linear_level:
                    if saturation_level:
                        if saturation_level > non_linear_level:
                            log.fullinfo(
                                f'Flagging non-linear pixels in {ad.filename} '
                                f'extension {ext.id} above level {non_linear_level:.2f}'
                            )
                            ext.mask |= np.where(
                                (ext.data >= non_linear_level)
                                & (ext.data < saturation_level),
                                DQ.non_linear,
                                0,
                            ).astype(DQ.datatype)
                            # Readout modes of IR detectors can result in
                            # saturated pixels having values below the
                            # saturation level. Flag those. Assume we have an
                            # IR detector here because both non-linear and
                            # saturation levels are defined and nonlin<sat
                            regions, nregions = measurements.label(
                                ext.data < non_linear_level
                            )
                            # In all my tests, region 1 has been the majority
                            # of the image; however, I cannot guarantee that
                            # this is always the case and therefore we should
                            # check the size of each region
                            region_sizes = measurements.labeled_comprehension(
                                ext.data,
                                regions,
                                np.arange(1, nregions + 1),
                                len,
                                int,
                                0,
                            )
                            # First, assume all regions are saturated, and
                            # remove any very large ones. This is much faster
                            # than progressively adding each region to DQ
                            hidden_saturation_array = np.where(
                                regions > 0, 4, 0
                            ).astype(DQ.datatype)
                            for region in range(1, nregions + 1):
                                # Limit of 10000 pixels for a hole is a bit arbitrary
                                if region_sizes[region - 1] > 10000:
                                    hidden_saturation_array[regions == region] = 0
                            ext.mask |= hidden_saturation_array

                        elif saturation_level < non_linear_level:
                            log.warning(
                                '%s extension %s has saturation level\
                                        less than non-linear level',
                                ad.filename,
                                ext.id,
                            )
                        else:
                            log.fullinfo(
                                'Saturation and non-linear levels are the same for\
                                         %s:%s. Only flagging saturated pixels',
                                ad.filename,
                                ext.id,
                            )
                    else:
                        log.fullinfo(
                            'Flagging non-linear pixels in %s:%s above level %s',
                            ad.filename,
                            ext.id,
                            non_linear_level,
                        )
                        ext.mask |= np.where(
                            ext.data >= non_linear_level, DQ.non_linear, 0
                        ).astype(DQ.datatype)
            if static and static.filename:
                add_provenance(
                    ad, static.filename, md5sum(static.path) or '', self.myself()
                )
            if user and user.filename:
                add_provenance(
                    ad, user.filename, md5sum(user.path) or '', self.myself()
                )

        # Handle latency if reqested
        if params.get('latency', False):
            try:
                adinputs = self.addLatencyToDQ(
                    adinputs, time=params['time'], non_linear=params['non_linear']
                )
            except AttributeError:
                log.warning(
                    'addLatencyToDQ() not defined in primitivesClass %s',
                    self.__class__.__name__,
                )

        # Add the illumination mask if requested
        if params['add_illum_mask']:
            adinputs = self.addIllumMaskToDQ(
                adinputs, **self._inherit_params(params, 'addIllumMaskToDQ')
            )

        # Timestamp and update filenames
        for ad in adinputs:
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=sfx, strip=True)

        return adinputs

    def validateData(self, adinputs=None, **params):
        """
        Validate MAROON-X data while ignoring invalid WCS exceptions.

        Parameters
        ----------
        adinputs : list of AstroData
            List of unchecked AstroData objects
        suffix : str
            suffix to be added to output files
        require_wcs : bool
            do all extensions have to have a defined WCS?

        Returns
        -------
        list of AstroData
            List of checked AstroData objects
        """
        try:
            super().validateData(adinputs, **params)
        except ValueError as e:
            if 'valid WCS' not in str(e):
                raise
        return adinputs

    def standardizeWCS(self, adinputs=None, **params):
        """MAROONXDR version of standarizeWCS to skip WCS processing."""
        log = self.log
        log.stdinfo('Skipping standarizeWCS() primitive.')
        return adinputs

    def checkArm(self, adinputs=None, **params):
        """
        Check arm consistency across all MAROON-X frames.

        Verify that all frames have the same camera arm (BLUE or RED) based on
        data tags. The first file sets the expected value. Currently assumes
        1 astrodata object comes from 1 single-extension FITS. Need to update
        if/when original FITS are MEF.

        Parameters
        ----------
        adinputs - list of un-checked MX frames

        Returns
        -------
        adoutputs - set of list that passes test,  always at least first frame
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))

        # find first object's MX-camera
        arm_set = (
            'BLUE'
            if 'BLUE' in adinputs[0].tags
            else 'RED'
            if 'RED' in adinputs[0].tags
            else 'UNDEFINED'
        )
        if arm_set == 'UNDEFINED':
            msg = f'{adinputs[0].filename} has no defined camera arm'
            log.error(msg)
            raise OSError(msg)
        adoutputs = []
        if len(adinputs) == 1:
            log.warning('Only one file passed to checkArm')

        # include other objects in list with same tag
        # Warn user and toss frame if not all frames are taken with the same arm
        for ad in adinputs:
            log.stdinfo(f"{ad.filename}: checking arm consistency")
            if arm_set not in ad.tags:
                log.warning(
                    'Not all frames taken with the same camera arm, '
                    'restricting set to first arm used in list'
                )
                warning_msg = f'Tossing frame: {ad.filename}'
                log.warning(warning_msg)
            else:
                ad.update_filename(suffix=params['suffix'], strip=True)
                adoutputs.append(ad)
        return adoutputs

    def checkMaster(self, adinputs=None, **params):
        """
        Check that MX frames are processed master frames.

        Parameters
        ----------
        adinputs - list of un-checked MX frames

        Returns
        -------
        adoutputs - set of list that passes test,  always at least first frame
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))

        # find first object's MX-camera
        master_type = (
            'DARK'
            if 'DARK' in adinputs[0].tags
            else 'FLAT'
            if 'FLAT' in adinputs[0].tags
            else 'UNDEFINED'
        )
        if master_type == 'UNDEFINED':
            msg = f'{adinputs[0].filename} unknown master type'
            log.error(msg)
            raise ValueError(msg)

        # these tags should define a master frame unique set of tags
        required_tags = {master_type, 'PROCESSED'}

        # include other objects in list with required tags
        # Warn user and toss frame if not all frames comply
        adoutputs = []
        for ad in adinputs:
            log.stdinfo(f"{ad.filename}: checking if master frame")
            if not required_tags.issubset(ad.tags):
                warning_msg = f'Tossing non-master {master_type} frame: {ad.filename}.'
                log.warning(warning_msg)
                continue

            ad.update_filename(suffix=params['suffix'], strip=False)
            adoutputs.append(ad)
        return adoutputs

    def correctImageOrientation(self, adinputs=None, **params):
        """
        Correct image orientation to proper echelle format.

        Flip SCI if needed so that left lower corner is bluest wavelength,
        upper right corner is reddest wavelength. Resulting echelle orders
        go from left to right. MX blue frames start with incorrect orientation
        for reduction. This primitive must be called after DQ is established
        and before any image arithmetic is performed.

        Parameters
        ----------
        adinputs - list of un-checked MX objects

        Returns
        -------
        adoutputs - same list as inputs, with correct orientation to SCI
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        timestamp_key = self.timestamp_keys[self.myself()]

        adoutputs = []
        for ad in adinputs:
            adout = deepcopy(ad)
            # check for tags inherited from original fits,
            # perform 2-axis flip as needed

            # Check if the frame was from the blue arm by looking at the tags
            # If it is, then flip the image
            # TODO: Change this to just look at the arm tag
            if 'BLUE' in ad.tags:
                log.stdinfo(f'{ad.filename} set as blue, orientation flipped')
                adout[0].data = np.fliplr(np.flipud(ad[0].data))
                try:
                    adout[0].mask = np.fliplr(np.flipud(ad[0].mask))
                except Exception:
                    log.warning(
                        'DQ not found for %s while orienting image', ad.filename
                    )

            # If it is not from the blue arm, then check if it is from the
            # red arm by looking at the image orientation. Do not flip the image
            elif 'RED' in ad.tags:
                log.stdinfo(f'{ad.filename} set as red, orientation unchanged')

            # In any other case, something has gone wrong- return an error
            else:
                msg = f'{ad.filename} has no defined orientation'
                log.error(msg)
                raise OSError(msg)
            adout.update_filename(suffix=params['suffix'], strip=True)
            adoutputs.append(adout)

        gt.mark_history(adoutputs, primname=self.myself(), keyword=timestamp_key)

        return adoutputs

    def addVAR(self, adinputs=None, **params):
        """
        Add variance extension to MAROON-X frames.

        Calculate the variance based on the read noise for the chip and the
        poisson noise (the variance in this case is just the number of photons
        for each pixel). The variance is then stored as a FITS extension for
        each file.

        Parameters
        ----------
        adinputs : list of AstroData
            List of MX objects without variance extensions
        suffix : str
            Suffix to be added to output files
        read_noise : bool
            Whether to include read noise in variance calculations
        poisson_noise : bool
            Whether to include poisson noise in variance calculations

        Returns
        -------
        list of AstroData
            List of MX objects with variance extensions
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        timestamp_key = self.timestamp_keys[self.myself()]

        read_noise = params['read_noise']
        poisson_noise = params['poisson_noise']

        adoutputs = []
        for ad in adinputs:
            log.stdinfo(f'Adding variance extension for {ad.filename}')

            npix_x, npix_y = ad[0].data.shape
            var = np.zeros((npix_x, npix_y), dtype=np.float32)

            if read_noise:
                # Check if the frame was from the blue arm by looking at the tags.
                # Use the read noise from first array Region
                read_noise_var = ad.read_noise()[0][0]
                var += read_noise_var

                log.stdinfo(
                    f'{ad.filename} read noise variance contribution is '
                    f'{read_noise_var}'
                )

            if poisson_noise:
                # The variance due to poisson noise is just the number of
                # photons for each pixel. Add to the read noise variance.
                var += ad[0].data

            ad[0].variance = var
            ad.update_filename(suffix=params['suffix'], strip=True)
            adoutputs.append(ad)
        gt.mark_history(adoutputs, primname=self.myself(), keyword=timestamp_key)
        return adoutputs

    def checkND(self, adinputs=None, **params):
        """
        Check ND filter consistency across all frames.

        Verify that the ND filter on the sim cal fiber is consistent through
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
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))

        # Get the simcal ND filter setting for the first file
        check_val = round(adinputs[0].filter_orientation()['ND'], 2)
        adoutputs = []

        # In case we have multiple files, check that they all have the same
        # ND filter setting
        if len(adinputs) > 1:
            for ad in adinputs:
                log.stdinfo(f"{ad.filename}: checking ND filter consistency")
                if check_val != round(ad.filter_orientation()['ND'], 2):
                    log.warning(
                        'Not all frames have the same simcal ND filter '
                        'setting, restricting set to first seen'
                    )
                else:
                    ad.update_filename(suffix=params['suffix'], strip=True)
                    adoutputs.append(ad)

            # If we only have one file with the correct filter setting, return an error
            if len(adoutputs) == 1:
                log.error(
                    'Only first frame found, of given, with its '
                    'simcal ND filter setting'
                )
                raise OSError

        # If we only have one file in total, return a warning
        #
        else:
            log.warning('Only one file passed to checkND')
            return adinputs
        return adoutputs

    def subtractOverscan(self, adinputs=None, **params):
        """
        Subtract overscan level from the image.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            log.stdinfo(f"{ad.filename}: subtracting overscan")
            # TODO: check if ovescan has already been subtracted

            # get the overscan section used for bias subtraction
            # and the array section where the correction is applied
            # MX has only one extension (red or blue), use index 0
            osec = ad.subtract_overscan_section()[0]
            asec = ad.array_subtract_overscan_section()[0]

            if len(osec) != len(asec):
                msg = 'Overscan and array sections for bias subtraction do not match.'
                log.warning(msg)
                raise ValueError(msg)

            # get data as float32
            ad[0].data = ad[0].data.astype(np.float32)
            data = ad[0].data

            # Subtract the mean overscan level from the array region
            for osec_, asec_ in zip(osec, asec):
                oslice = osec_.asslice()
                aslice = asec_.asslice()
                data[aslice] -= np.mean(data[oslice])  # TODO: use nanmean

            # Timestamp, and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params['suffix'], strip=True)

        return adinputs

    def stackFramesMXCal(self, adinputs=None, **params):
        """
        Stack MAROON-X calibration frames with etalon flux scaling.

        MX-specific version of stackFrames for calibration frames - changes
        scaling to average full frame mean to purposely scale by etalon flux
        and its drift between calibration exposures. This function should
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

        def flatten_item(*args):
            return (
                el
                for item in args
                for el in (
                    flatten_item(*item) if isinstance(item, (list, tuple)) else (item,)
                )
            )

        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        timestamp_key = self.timestamp_keys['stackFrames']
        sfx = params['suffix']
        memory = params['memory']
        if memory is not None:
            memory = int(memory * 1000000000)

        zero = params['zero']
        scale = params['scale']
        apply_dq = params['apply_dq']
        separate_ext = params['separate_ext']
        statsec = params['statsec']
        reject_method = params['reject_method']
        save_rejection_map = params['save_rejection_map']

        if statsec:
            statsec = tuple(
                [
                    slice(int(start) - 1, int(end))
                    for x in reversed(statsec.strip('[]').split(','))
                    for start, end in [x.split(':')]
                ]
            )

        # Check that the input AstroData objects are compatible
        if len(adinputs) <= 1:
            log.stdinfo(
                'No stacking will be performed, since at least two '
                'input AstroData objects are required for stackFrames'
            )
            return adinputs

        if (
            reject_method == 'minmax'
            and self.mode == 'qa'
            and params['nlow'] + params['nhigh'] >= len(adinputs)
        ):
            log.warning('Trying to reject too many images. Setting nlow=nhigh=0.')
            params['nlow'] = 0
            params['nhigh'] = 0

        if len({len(ad) for ad in adinputs}) > 1:
            msg = 'Not all inputs have the same number of extensions'
            log.debug(msg)
            raise OSError(msg)
        if len({ext.nddata.shape for ad in adinputs for ext in ad}) > 1:
            msg = 'Not all inputs images have the same shape'
            log.debug(msg)
            raise OSError(msg)

        num_img = len(adinputs)
        num_ext = len(adinputs[0])
        zero_offsets = np.zeros((num_ext, num_img), dtype=np.float32)
        scale_factors = np.ones_like(zero_offsets)

        # Try to determine how much memory we're going to need to stack and
        # whether it's necessary to flush pixel data to disk first
        # Also determine kernel size from offered memory and bytes per pixel
        num_bytes_per_ext = []
        for ext in adinputs[0]:
            num_bytes = 0
            # Count _data twice to handle temporary arrays
            num_bytes += 2 * ext.data.dtype.itemsize
            num_bytes += 2  # mask always created
            num_bytes_per_ext.append(num_bytes * np.prod(ext.shape))

        if memory is not None and (num_img * max(num_bytes_per_ext) > memory):
            adinputs = self.flushPixels(adinputs)

        # Compute the scale and offset values by accessing the memmapped data
        # so we can pass those to the stacking function
        # TODO: Should probably be done better to consider only the overlap
        # regions between frames
        if scale or zero:
            levels = np.empty((num_img, num_ext), dtype=np.float32)
            for i, ad in enumerate(adinputs):
                for index in range(num_ext):
                    nddata = (
                        ad[index].nddata.window[:]
                        if statsec is None
                        else ad[index].nddata.window[statsec]
                    )
                    scale_mask = (
                        ad[index].nddata.window[:].mask
                        if statsec is None
                        else ad[index].nddata.window[statsec].mask
                    ) == DQ.good
                    # MX specific changed line
                    # uses entire good pixel frame, purposely
                    # including etalon flux, to calculate level as sum
                    levels[i, index] = np.nansum(nddata.data[scale_mask])

            if scale and zero:
                log.warning('Both scale and zero are set. Setting scale=False.')
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
                    scale_factors = np.tile(
                        target / np.mean(levels, axis=1), num_ext
                    ).reshape(num_ext, num_img)
                else:
                    zero_offsets = np.tile(
                        target - np.mean(levels, axis=1), num_ext
                    ).reshape(num_ext, num_img)

            # Check for negative, infinite or undefined scale factors
            if scale and np.min(scale_factors) < 0:
                log.warning('Some scale factors are negative. Not scaling.')
                scale_factors = np.ones_like(scale_factors)
                scale = False
            if scale and np.any(np.isinf(scale_factors)):
                log.warning('Some scale factors are infinite. Not scaling.')
                scale_factors = np.ones_like(scale_factors)
                scale = False
            if scale and np.any(np.isnan(scale_factors)):
                log.warning('Some scale factors are undefined. Not scaling.')
                scale_factors = np.ones_like(scale_factors)
                scale = False

        if reject_method == 'varclip' and any(
            ext.variance is None for ad in adinputs for ext in ad
        ):
            log.warning(
                "Rejection method 'varclip' has been chosen but some"
                " extensions have no variance. 'sigclip' will be used"
                ' instead.'
            )
            reject_method = 'sigclip'

        log.stdinfo(
            f"Combining {num_img} inputs with {params['operation']}\
                     and {reject_method} rejection"
        )

        stack_function = NDStacker(
            combine=params['operation'], reject=reject_method, log=self.log, **params
        )

        # NDStacker uses DQ if it exists; if we don't want that, delete the DQs!
        if not apply_dq:
            for ad in adinputs:
                for ext in ad:
                    ext.mask = None  # delete mask

        ad_out = astrodata.create(adinputs[0].phu)
        for index, (ext, sfactors, zfactors) in enumerate(
            zip(adinputs[0], scale_factors, zero_offsets)
        ):
            status = (
                f'Combining extension {ext.id}.' if num_ext > 1 else 'Combining images.'
            )
            if scale:
                status += ' Applying scale factors.'
                numbers = sfactors
            elif zero:
                status += ' Applying offsets.'
                numbers = zfactors
            log.stdinfo(status)
            if (scale or zero) and (index == 0 or separate_ext):
                for ad, value in zip(adinputs, numbers):
                    # need one digit beyond 10.3f to see differences
                    log.stdinfo(f'{ad.filename:40s}{value:10.4f}')

            shape = adinputs[0][index].nddata.shape
            if memory is None:
                kernel = shape
            else:
                # Chop the image horizontally into equal-sized chunks to process
                # This uses the minimum number of steps and uses minimum memory
                # per step.
                oversubscription = (num_bytes_per_ext[index] * num_img) // memory + 1
                kernel = (
                    (shape[0] + oversubscription - 1) // oversubscription,
                ) + shape[1:]

            with_mask = apply_dq and not any(
                ad[index].nddata.window[:].mask is None for ad in adinputs
            )
            result = windowedOp(
                stack_function,
                [ad[index].nddata for ad in adinputs],
                scale=sfactors,
                zero=zfactors,
                kernel=kernel,
                dtype=np.float32,
                with_uncertainty=True,
                with_mask=with_mask,
                save_rejection_map=save_rejection_map,
            )
            ad_out.append(result)
            log.stdinfo('')

        # Set AIRMASS to be the mean of the input values
        try:
            airmass_kw = ad_out._keyword_for('airmass')
            mean_airmass = np.mean([ad.airmass() for ad in adinputs])
        except Exception:  # generic implementation failure (probably non-Gemini)
            pass
        else:
            ad_out.phu.set(airmass_kw, mean_airmass, 'Mean airmass for the exposure')

        # Add suffix to datalabel to distinguish from the reference frame
        if sfx[0] == '_':
            extension = sfx.replace('_', '-', 1).upper()
        else:
            extension = '-' + sfx.upper()
        ad_out.phu.set(
            'DATALAB',
            f'{ad_out.data_label()}{extension}',
            self.keyword_comments['DATALAB'],
        )

        # Add other keywords to the Fits Header Unit about the stacking inputs
        ad_out.orig_filename = ad_out.phu.get('ORIGNAME')
        ad_out.phu.set(
            'NCOMBINE %d %s', len(adinputs), self.keyword_comments['NCOMBINE']
        )
        for i, ad in enumerate(adinputs, start=1):
            ad_out.phu.set(f'IMCMB{i:03d}', ad.phu.get('ORIGNAME', ad.filename))

        # Timestamp and update filename and prepare to return single output
        gt.mark_history(ad_out, primname=self.myself(), keyword=timestamp_key)
        ad_out.update_filename(suffix=sfx, strip=True)

        return [ad_out]

    def stackDarks(self, adinputs=None, **params):
        """
        Stack MAROON-X dark frames with etalon intensity scaling.

        MX-specific version of stack darks allowing scaling for etalon
        intensity drift that is in MX 'darks'.

        Parameters
        ----------
        adinputs : list of AstroData
            Input frames to be combined
        suffix : str
            Suffix to be added to output files
        scale_mode : str
            Scaling mode for the input frames.
            Options are 'mean_frame' or 'first_frame'.
        lsigma : float
            Lower sigma clipping threshold for the rejection method
        hsigma : float
            Upper sigma clipping threshold for the rejection method
        max_iters : int
            Maximum number of iterations for the rejection method
        reject_method : str
            Rejection method to be used. Currently only 'sigclip' is supported.

        Returns
        -------
        list of AstroData
            Combined output frame
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        timestamp_key = self.timestamp_keys['stackFrames']

        # Check that all inputs are DARKs and have same exposure time
        if not all('DARK' in dark.tags for dark in adinputs):
            msg = 'Not all inputs have DARK tag'
            log.debug(msg)
            raise ValueError(msg)

        if not all(
            dark.exposure_time() == adinputs[0].exposure_time() for dark in adinputs[1:]
        ):
            msg = 'Darks are not of equal exposure time'
            log.debug(msg)
            raise ValueError(msg)

        if len(adinputs) <= 1:
            log.stdinfo(
                'No stacking will be performed, since at least two '
                'input AstroData objects are required for stackFrames'
            )
            return adinputs

        sfx = params['suffix']
        scale_mode = params['scale_mode']
        lsigma = params['lsigma']
        hsigma = params['hsigma']
        max_iters = params['max_iters']
        reject_method = params['reject_method']

        # Create output frame structure
        ad_out = deepcopy(adinputs[0])

        # STACKING BEGINS ================================================
        # Stack frames into a data cube
        data_cube = np.dstack([ad[0].data.astype(np.float32) for ad in adinputs])

        # Scale frames
        data_cube = _scaleCube(data_cube, scale_mode=scale_mode)

        # Combine frames using the specified operation
        if reject_method != 'sigclip':
            msg = f'Unsupported rejection method: {reject_method}'
            log.debug(msg)
            raise ValueError(msg)

        data_mean, data_median, data_stddev = sigma_clipped_stats(
            data_cube,
            axis=2,
            sigma_lower=lsigma,
            sigma_upper=hsigma,
            maxiters=max_iters,
        )

        # Further bad pixel filtering
        # this logic is inherited from the legacy code
        data_min = np.min(data_cube, axis=2)
        diff = data_mean - data_min
        bad_pixels = np.where(
            (diff > 2 * np.nanmedian(data_mean))
            & (diff > np.abs(0.3 * data_min))
            & (data_mean > 10)
            & (diff > 10)
        )
        n_bad = np.shape(bad_pixels)[1]
        if n_bad > 0:
            data_mean[bad_pixels] = data_min[bad_pixels]

        # Set the output data
        # Dragons enforces float32.
        ad_out[0].data = data_mean.astype(np.float32)
        # ================================================================

        # Add suffix to datalabel to distinguish from the reference frame
        if sfx[0] == '_':
            extension = sfx.replace('_', '-', 1).upper()
        else:
            extension = '-' + sfx.upper()
        ad_out.phu.set(
            'DATALAB',
            f'{ad_out.data_label()}{extension}',
            self.keyword_comments['DATALAB'],
        )

        # Add other keywords to the PHU about the stacking inputs
        ad_out.orig_filename = ad_out.phu.get('ORIGNAME')
        ad_out.phu.set('NCOMBINE', len(adinputs), self.keyword_comments['NCOMBINE'])
        for i, ad in enumerate(adinputs, start=1):
            ad_out.phu.set(f'IMCMB{i:03d}', ad.phu.get('ORIGNAME', ad.filename))

        # Timestamp and update filename and prepare to return single output
        gt.mark_history(ad_out, primname=self.myself(), keyword=timestamp_key)
        ad_out.update_filename(suffix=sfx, strip=True)

        return [ad_out]

    def stackDarksOld(self, adinputs=None, **params):
        """
        Stack MAROON-X dark frames using old method.

        MX-specific version of stack darks allowing scaling for etalon
        intensity drift that is in MX 'darks'.
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))

        # Check that all inputs are DARKs and have same exposure time
        if not all('DARK' in dark.tags for dark in adinputs):
            msg = 'Not all inputs have DARK tag'
            log.debug(msg)
            raise ValueError(msg)

        if not all(
            dark.exposure_time() == adinputs[0].exposure_time() for dark in adinputs[1:]
        ):
            msg = 'Darks are not of equal exposure time'
            log.debug(msg)
            raise ValueError(msg)

        # MX specific-changed lines start
        # MX 'dark' frames have flux in them, need to scale.
        # Also utilizes special stackFramesMXCal scaling.
        stack_params = self._inherit_params(params, 'stackFramesMXCal')
        stack_params.update({'zero': False})
        adinputs = self.stackFramesMXCal(adinputs, **params)
        # MX specific-changed lines end

        return adinputs

    def stackFlatsOld(self, adinputs=None, **params):
        """
        Stack MAROON-X flat frames using old method.

        MaroonX-specific version of stack flats to call correct stackframes
        for the flats.
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        # MX specific-changed line
        # Utilizes special stackFramesMXCal scaling.
        stack_params = self._inherit_params(params, 'stackFramesMXCal')
        stack_params.update({'zero': False})
        adinputs = self.stackFramesMXCal(adinputs, **stack_params)
        return adinputs

    def stackFlats(self, adinputs=None, **params):
        """
        Stack MAROON-X flat field frames.

        A simplified DRAGONS primitive that reproduces the legacy make_master_flat.py.

        Parameters
        ----------
        adinputs : list of AstroData
            Input frames to be combined
        operation : str
            Combine method (default: 'mean')
        reject_method : str
            Rejection method (default: 'sigclip')
        scale : bool
            Whether to scale images (default: True)

        Returns
        -------
        list of AstroData
            Combined output frame
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        timestamp_key = self.timestamp_keys['stackFrames']

        sfx = params['suffix']
        scale_mode = params['scale_mode']
        lsigma = params['lsigma']
        hsigma = params['hsigma']
        max_iters = params['max_iters']
        reject_method = params['reject_method']

        if len(adinputs) <= 1:
            log.stdinfo(
                'No stacking will be performed, since at least two '
                'input AstroData objects are required for stackFrames'
            )
            return adinputs

        # Check that all inputs are FLAT
        if not all('FLAT' in ad.tags for ad in adinputs):
            msg = 'Not all inputs have FLAT tag'
            log.debug(msg)
            raise ValueError(msg)

        # STACKING BEGINS ================================================
        # Create output frame structure
        ad_out = deepcopy(adinputs[0])

        # Stack frames into a data cube
        data_cube = np.dstack([ad[0].data.astype(np.float32) for ad in adinputs])

        # Scale frames
        data_cube = _scaleCube(data_cube, scale_mode=scale_mode)

        # Combine frames using the specified operation
        if reject_method != 'sigclip':
            msg = f'Unsupported rejection method: {reject_method}'
            log.debug(msg)
            raise ValueError(msg)

        data_mean, data_median, data_stddev = sigma_clipped_stats(
            data_cube,
            axis=2,
            sigma_lower=lsigma,
            sigma_upper=hsigma,
            maxiters=max_iters,
        )

        # Further bad pixel filtering
        # this logic is inherited from the legacy code
        data_mean[data_mean < 0] = 0

        # sigma_clipped_stats promotes float32 input to float64.
        # Legacy keeps the float64 result (BITPIX=-64) at this stage;
        # Dragons enforces float32.
        ad_out[0].data = data_mean.astype(np.float32)
        # ================================================================

        # Add suffix to datalabel to distinguish from the reference frame
        if sfx[0] == '_':
            extension = sfx.replace('_', '-', 1).upper()
        else:
            extension = '-' + sfx.upper()
        ad_out.phu.set(
            'DATALAB',
            f'{ad_out.data_label()}{extension}',
            self.keyword_comments['DATALAB'],
        )

        # Add other keywords to the PHU about the stacking inputs
        ad_out.orig_filename = ad_out.phu.get('ORIGNAME')
        ad_out.phu.set('NCOMBINE', len(adinputs), self.keyword_comments['NCOMBINE'])
        for i, ad in enumerate(adinputs, start=1):
            ad_out.phu.set(f'IMCMB{i:03d}', ad.phu.get('ORIGNAME', ad.filename))

        # Timestamp and update filename and prepare to return single output
        gt.mark_history(ad_out, primname=self.myself(), keyword=timestamp_key)
        ad_out.update_filename(suffix=sfx, strip=True)

        return [ad_out]

    def findStripes(self, adinputs=None, **params):
        """
        Locate and fit stripes in flat field spectrum.

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
        ----------
        adinputs : list of AstroData
            Single MX astrodata object, is either a DFFFD flat,
            FDDDF flat, or combined FFFFF flat
        suffix : str
            Suffix to be added to output files
        deg_polynomial : int
            Degree of the polynomial fit
        med_filter : int
            Median filter parameter
        gauss_filter_sigma : float
            Sigma of the gaussian filter used to smooth the image
        min_peak : float
            Minimum relative peak height

        Returns
        -------
        list of AstroData
            Single MX astrodata object with STRIPES_LOC extension.
            This extension temporarily holds the fits-unsavable fiber information
            before it is utilized and then removed.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        deg_polynomial = params['deg_polynomial']
        med_filter = params['med_filter']
        gauss_filter_sigma = params['gauss_filter_sigma']
        min_peak = params['min_peak']

        for ad in adinputs:
            log.stdinfo(f"{ad.filename}: finding stripes")
            npix_y, npix_x = ad[0].data.shape
            # first smooth frame to improve algorithm stability
            filtered_ad = median_filter(ad[0].data, med_filter)
            filtered_ad = gaussian_filter(filtered_ad, gauss_filter_sigma)

            # find peaks in center column
            data = filtered_ad[:, int(npix_x / 2)]
            peaks = (
                np.r_[True, data[1:] >= data[:-1]] & np.r_[data[:-1] > data[1:], True]
            )

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
            log.fullinfo(f'Fit polynomial of order {deg_polynomial} ' f'to each stripe')
            x_pixels = np.arange(npix_x)
            polynomials = [
                np.poly1d(np.polyfit(x_pixels, o, deg_polynomial)) for o in orders
            ]
            # polynomials is FITs-unsaveable dict of dicts,
            # will use and then remove before storing
            ad[0].STRIPES_LOC = polynomials
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params['suffix'], strip=False)
        return adinputs

    def identifyStripes(self, adinputs=None, **params):
        """
        Identify stripes by assigning order and fiber numbers.

        Assign proper order and fiber number to each stripe, including
        correction for the possibility that the spectra have shifted up/down
        in the cross-dispersion direction since the reference was made.
        Requires the findStripes primitive to be run prior during recipe so
        the stripes are located in the input, i.e. STRIPES_LOC extension exists.

        Parameters
        ----------
        adinputs : list of AstroData
            Single MX astrodata object, is either a DFFFD flat,
            FDDDF flat, or combined FFFFF flat with STRIPES_LOC extension
        suffix : str
            Suffix to be added to output files
        positions_dir : str
            Lookup fits location of nominal y positions and
            fiber/order labels. Shape is Nx3, columns are
            [fibers, orders, y_positions], nominally found in lookups/SID
        selected_fibers : list of int
            List of fiber numbers illuminated in the flat, if None assumes all
            can work if not given on partially illuminated frame, but best
            practice is to explicitly identify on function call

        Returns
        -------
        list of AstroData
            Single MX astrodata object with STRIPES_ID extension.
            This extension temporarily holds the fits-unsavable fiber information
            before it is utilized and then removed. A new extension REMOVED_STRIPES
            also saves the polynomial info for every stripe that is not identified
            from the original set inherited from findStripes
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        positions_dir = params.get('positions_dir', None)
        selected_fibers = params.get('selected_fibers', None)

        for ad in adinputs:
            log.stdinfo(f"{ad.filename}: identifying stripes")
            # first obtain reference positions (position_dir)
            if positions_dir is None:
                positions_dir = maroonx_utils.get_sid_filename(ad)
                log.fullinfo(f"Using SID file: {positions_dir}")
            with fits.open(positions_dir, 'readonly') as lookup_frame:
                positions = lookup_frame[1].data
            _, npix_x = ad.data[0].shape
            p_id = {}

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
                    closest_stripe_idx = np.argmin(np.abs(y_positions + shift - y))
                    distances.append(
                        np.abs(y_positions[closest_stripe_idx] + shift - y)
                    )
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
            for _, poly in enumerate(ad[0].STRIPES_LOC):  # for each stripe
                # Observed y position of stripe is given by the stored
                # polynomial and the center column of frame
                observed_y = np.poly1d(poly)(npix_x / 2)
                # Find the closest stripe in the reference frame,
                # using the SID file's y positions
                closest_stripe_idx = np.argmin(
                    np.abs(y_positions + shift_calculated - observed_y)
                )

                # If the stripe is within 7 pixels of the reference, assign it
                if (
                    np.abs(
                        y_positions[closest_stripe_idx] + shift_calculated - observed_y
                    )
                    < 7
                ):
                    if used[closest_stripe_idx] == 0:
                        used[closest_stripe_idx] = 1
                        fiber = fibers[closest_stripe_idx]
                        order = orders[closest_stripe_idx]
                        log.fullinfo('fiber %s, order %s found', fiber, order)

                        fiber = f'fiber_{fiber}'
                        order = f'{order}'

                        if fiber in p_id:
                            p_id[fiber].update({order: poly.coefficients})
                        else:
                            p_id[fiber] = {order: poly.coefficients}
                    # If the stripe is not in the dictionary, warn the user
                    else:
                        log.warning(
                            'Stripe at %s could not be identified unambiguously',
                            observed_y,
                        )
                        unidentified_stripes.append(poly.coefficients)
                # If stripe not within 7 pixels, warn and add to unidentified
                else:
                    log.warning('Stripe at %s could not be identified', observed_y)
                    unidentified_stripes.append(poly.coefficients)
            # p_id is FITs-unsavable dict of dicts,
            # will use and reformat before storing
            ad[0].STRIPES_ID = p_id
            # REMOVED_STRIPES is FITS-savable np.array of unidentified
            # stripes with no meta-info loss as stripe is unidentified
            ad[0].REMOVED_STRIPES = np.array(unidentified_stripes)
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params['suffix'], strip=False)
        return adinputs

    def defineFlatStripes(self, adinputs=None, **params):
        """
        Define flat field stripe locations for extraction and straylight removal.

        Save fiber location info based on flat field info for stray light
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
        ----------
        adinputs : list of AstroData
            Single MX astrodata object, is either a DFFFD, FDDDF flat,
            or combined FFFFF flat
        suffix : str
            Suffix to be added to output files
        slit_height : int
            Half pixel height of box in each dimension to
            perform box extraction with
        extract : bool
            If True, will write STRIPES_ID in fits-acceptable
            format. Utilized in combined, all fiber illuminated FFFFF_flat

        Returns
        -------
        list of AstroData
            Single MX astrodata object with INDEX_FIBER, INDEX_ORDER
            extensions and possibly STRIPES_ID and STRIPES_FIBERS extensions
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        slit_height = params['slit_height']
        extract = params['extract']

        for ad in adinputs:
            log.stdinfo(f"{ad.filename}: defining flat stripes")
            p_id = ad[0].STRIPES_ID  # Create dictionary of stripes
            img = ad[0].data
            npix_y, npix_x = img.shape
            x_pixels = np.arange(npix_x)
            index_fiber = np.zeros_like(img, dtype=np.int8)  # 2D pixel map
            index_order = np.zeros_like(img, dtype=np.int8)
            slit_indices_y = (
                np.arange(-slit_height, slit_height)
                .repeat(npix_x)
                .reshape((2 * slit_height, npix_x))
            )
            slit_indices_x = np.tile(np.arange(npix_x), 2 * slit_height).reshape(
                (2 * slit_height, npix_x)
            )
            for fiber_iter in p_id.keys():
                for order_iter, poly in p_id[fiber_iter].items():
                    y = np.poly1d(poly)(x_pixels)
                    indices = np.rint(slit_indices_y + y).astype(int)
                    valid_indices = np.logical_and(indices < npix_y, indices > 0)
                    final_fiber = fiber_iter
                    if isinstance(final_fiber, str):
                        final_fiber = int(''.join(filter(str.isdigit, final_fiber)))
                    index_fiber[
                        indices[valid_indices], slit_indices_x[valid_indices]
                    ] = final_fiber
                    final_order = order_iter
                    if isinstance(final_order, str):
                        final_order = int(''.join(filter(str.isdigit, final_order)))
                    index_order[
                        indices[valid_indices], slit_indices_x[valid_indices]
                    ] = final_order

            ad[0].INDEX_FIBER = index_fiber.astype(int)
            ad[0].INDEX_ORDER = index_order.astype(int)
            del ad[0].STRIPES_ID  # delete interim information in improper format
            del ad[0].STRIPES_LOC  # delete interim information

            if extract:
                fiber_tables = []
                for ifib in sorted(p_id.keys(), key=lambda x: x.lower()):
                    fiber_tables.append(Table.from_pandas(pd.DataFrame(p_id[ifib])))
                ad[0].STRIPES_ID = vstack(fiber_tables, metadata_conflicts='silent')
                ad[0].STRIPES_FIBERS = np.array(
                    [key[-1] for key in sorted(p_id.keys(), key=lambda x: x.lower())]
                ).astype(int)
                # actually extracting stripes using full info in DRAGONS,
                # aka sparse matrix realizations,
                # not utilized as csc not FITS compatible

            ad.update_filename(suffix=params['suffix'], strip=True)
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)
        return adinputs

    def removeStrayLight(self, adinputs=None, **params):
        """
        Remove stray light from full frame images.

        Remove stray light from full frame images for more accurate fiber flux
        accounting. Requires the defineStripes primitive to be run prior during
        recipe so INDEX_FIBER and INDEX_ORDER extensions exist to define pixel
        locations across frame within fiber traces to avoid when finding stray
        light.

        Parameters
        ----------
        adinputs : list of AstroData
            Single MX astrodata object, is either a DFFFD or FDDDF flat
            that has not previously had its stray light removed
        suffix : str
            Suffix to be added to output files
        box_size : int
            Pixel height and width of 'mesh_element' used in
            background identification sub-routine
        filter_size : int
            Pixel height and width of window to perform
            background identification sub-routine
        snapshot : bool
            Bool to save difference frame of removed stray light as
            extension STRAYLIGHT_DIFFERENCE
        report : bool
            Generate PDF diagnostic report showing straylight removal stages

        Returns
        -------
        list of AstroData
            Single MX astrodata object with stray light removed from SCI
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        if params.pop('legacy'):
            # Run the legacy patch version and exit
            log.fullinfo('Running Legacy Patch removeStrayLight')
            adoutputs = self.removeStrayLight_legacyPatch(adinputs, **params)
            return adoutputs

        box_size = params['box_size']
        filter_size = params['filter_size']
        snapshot = params['snapshot']
        report = params['report']

        adoutputs = []
        for ad in adinputs:
            log.stdinfo(f"{ad.filename}: removing straylight")
            adint = deepcopy(ad)
            adint2 = deepcopy(ad)
            adout = deepcopy(ad)
            stripes = ad[0].INDEX_FIBER
            orders = ad[0].INDEX_ORDER

            # Mask all actual bad pixels from bpm
            adint[0].data[adint[0].mask != DQ.good] = np.nan
            max_val_masked = np.nanmax(adint[0].data)

            # Mask all the stripes so we only fit the background
            adint[0].data[stripes > 0] = np.nan

            # Check if we are on the blue chip
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

            # Else we must be on the red chip
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
            bkg = Background2D(
                adint[0].data,
                (box_size, box_size),
                filter_size=(filter_size, filter_size),
                sigma_clip=SigmaClip(sigma=4.0),
                bkg_estimator=MedianBackground(),
                exclude_percentile=95,
            )

            # in the blue, we may overshoot in the region where stripes are
            # missing or orders >121 a second round of background fitting
            # for negative and near-negative pixels should correct that

            adint2[0].data = ad[0].data - bkg.background
            median_negative = np.median(adint2[0].data[adint2[0].data < 0])
            log.fullinfo(
                f'Median sub-zero pixel value after 1st iteration: '
                f'{median_negative}'
            )
            adint2[0].data[adint2[0].data > 2] = np.nan
            adint2[0].data[adint2[0].data < -50] = np.nan
            bkg2 = Background2D(
                adint2[0].data,
                (box_size, box_size),
                filter_size=(filter_size, filter_size),
                sigma_clip=SigmaClip(sigma=4.0),
                bkg_estimator=MedianBackground(),
                exclude_percentile=95,
            )

            # Perform the second background fitting
            adout[0].data = ad[0].data - bkg.background - bkg2.background
            log.fullinfo(
                f'Truncated median sub-zero pixel value: '
                f'{np.median(adout[0].data[adout[0].data < 0])}'
            )

            # Generate PDF report if requested
            if report:
                # Capture result before clipping negative values
                result_before_clip = adout[0].data.copy()
                pdf_name = ad.filename.replace('.fits', '_backgroundfit.pdf')
                figs = plot_backgroundfit(
                    ad[0].data, adint[0].data, bkg, bkg2,
                    result_before_clip, max_val_masked,
                )
                with PdfPages(pdf_name) as pdf:
                    for fig in figs:
                        pdf.savefig(fig)
                        plt.close(fig)
                log.stdinfo(f"Straylight diagnostic report written to '{pdf_name}'")

            adout[0].data[adout[0].data < 0] = 0.01

            if snapshot:
                # used for unit testing, save straylight difference
                # but return original data so it can be compared in test
                adout[0].STRAYLIGHT_DIFFERENCE = adout[0].data - deepcopy(ad)[0].data
                adout[0].data = deepcopy(ad)[0].data
            adout.update_filename(suffix=params['suffix'], strip=True)
            adoutputs.append(adout)

        gt.mark_history(adoutputs, primname=self.myself(), keyword=timestamp_key)

        return adoutputs

    def removeStrayLight_legacyPatch(self, adinputs=None, **params):
        """
        Remove stray light from full frame images using legacy arrays.

        Remove stray light from full frame images for more accurate fiber flux
        accounting. Requires the defineStripes primitive to be run prior during
        recipe so INDEX_FIBER and INDEX_ORDER extensions exist to define pixel
        locations across frame within fiber traces to avoid when finding stray
        light.

        Parameters
        ----------
        adinputs : list of AstroData
            Single MX astrodata object, is either a DFFFD or FDDDF flat
            that has not previously had its stray light removed
        suffix : str
            Suffix to be added to output files
        box_size : int
            Pixel height and width of 'mesh_element' used in
            background identification sub-routine
        filter_size : int
            Pixel height and width of window to perform
            background identification sub-routine
        snapshot : bool
            Bool to save difference frame of removed stray light as
            extension STRAYLIGHT_DIFFERENCE

        Returns
        -------
        list of AstroData
            Single MX astrodata object with stray light removed from SCI
        """
        log = self.log
        log.debug(gt.log_message("primitive", "removeStrayLight", "starting"))
        timestamp_key = self.timestamp_keys['removeStrayLight']

        snapshot = params['snapshot']
        report = params['report']

        from pathlib import Path
        import os
        
        legacy_test = os.environ.get('MAROONX_LEGACY_TEST')
        if not legacy_test:
            raise RuntimeError(
                'MAROONX_LEGACY_TEST environment variable is not set.'
            )
        legacy_path = Path(legacy_test).parent / 'legacy_bkg_arrays'
        if not legacy_path.exists():
            raise FileNotFoundError(
                f'Legacy bkg arrays directory not found: {legacy_path}'
            )

        def _legacy_straylight_filename(base, legacy_path):
            parts = base.split('_')
            timestamp = parts[0]
            fiber_cfg = parts[1]
            arm = parts[2]
            seq = parts[3]

            if fiber_cfg in ('DFFFD', 'DDDDF'):
                short_ts = timestamp[:11]  # 'YYYYMMDDTHH', e.g. '20241114T18'
                name = f'{short_ts}_masterflat_{fiber_cfg}_{arm}_{seq}.npy'
            elif fiber_cfg == 'SOOOE':
                name = f'{timestamp}_{fiber_cfg}_{arm}_{seq}_backgroundfit_straylight.npy'
            else:
                raise ValueError(f'Unknown fiber configuration: {fiber_cfg}')

            return legacy_path / name

        adoutputs = []
        for ad in adinputs:
            log.stdinfo(f"{ad.filename}: removing straylight (legacy patch)")
            # Remove extension from outfile
            base = Path(ad.filename).stem.removesuffix('_overscanSubtracted')
            # data_dict = np.load(file_dict[base], allow_pickle=True).item()
            file_name = _legacy_straylight_filename(base, legacy_path)
            data_dict = np.load(file_name, allow_pickle=True).item()

            adout = deepcopy(ad)

            # Check we are loading the correct file.
            # Compare at float32: the .npy snapshots were recorded from
            # float32 data; the pipeline now keeps float64 until this point.
            np.testing.assert_allclose(
                data_dict['data'].astype(np.float32),
                ad[0].data.astype(np.float32),
            )

            # -------------------------------------------------------------------
            # Replace with legacy results.  Legacy backgroundfit.py reads
            # stacked flats as float32 (line 237), so the result is float32.
            # adout[0].data = data_dict['result'].astype(np.float32)
            
            # Flat frames go through a FITS roundtrip in legacy (BITPIX=-32), truncating to float32.
            # Science frames do not — legacy keeps float64 all the way to extraction.
            if 'FLAT' in ad.tags:
                adout[0].data = data_dict['result'].astype(np.float32)
            else:
                adout[0].data = data_dict['result']
            # -------------------------------------------------------------------
            if snapshot:
                adout[0].STRAYLIGHT_DIFFERENCE = adout[0].data - ad[0].data
                adout[0].data = ad[0].data

            # Generate PDF report if requested
            if report:
                # Capture result before clipping negative values
                result_before_clip = adout[0].data.copy()
                pdf_name = ad.filename.replace('.fits', '_backgroundfit_legacy.pdf')
                figs = plot_backgroundfit(
                    data_dict['data'], data_dict['data_masked'],
                    data_dict['bkg_background'], data_dict['bkg2_background'],
                    result_before_clip, data_dict['max_val_masked'],
                    bkg1_mesh_nmasked=data_dict['bkg_mesh_nmasked'],
                    box_size=data_dict['box_size'],
                )
                with PdfPages(pdf_name) as pdf:
                    for fig in figs:
                        pdf.savefig(fig)
                        plt.close(fig)
                log.stdinfo(f"Straylight diagnostic report written to '{pdf_name}'")

            adout.update_filename(suffix=params['suffix'], strip=True)
            adoutputs.append(adout)

        gt.mark_history(adoutputs, primname='removeStrayLight', keyword=timestamp_key)

        return adoutputs

    def separateFlatStreams(self, adinputs=None, **params):
        """
        Separate flat data into two streams based on fiber setup.

        This primitive divides the input flats into two categories:
        - 'DFFFD': stored in p.streams['DFFFD_flats']
        - 'FDDDF' or 'DDDDF': stored in p.streams['main']

        Parameters
        ----------
        adinputs : list of MX flats
        **params : dict of parameters

        Returns
        -------
        adoutputs : 'FDDDF' or 'DDDDF' list of AstroData objects

        Notes
        -----
        Modifies the instance's `streams` dictionary by adding a 'DFFFD_flats' key
        containing the list of DFFFD flat field AstroData objects.

        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))

        # Define flat setups to sort the streams
        DARK, FLAT = 'Dark', 'Flat lamp'
        FDDDF_setup = [FLAT, DARK, DARK, DARK, FLAT]
        DDDDF_setup = [DARK, DARK, DARK, DARK, FLAT]
        DFFFD_setup = [DARK, FLAT, FLAT, FLAT, DARK]

        # Initialize lists of AstroData objects to be added to the streams
        flat_fdddf_list = []
        flat_dfffd_list = []
        mislabeled = []

        for ad in adinputs:
            # Create list of FDDDF flats to go in the main stream
            if 'FLAT' in ad.tags and ad.fiber_setup() in [FDDDF_setup, DDDDF_setup]:
                flat_fdddf_list.append(ad)
                log.fullinfo(f'FDDDF Flat: {ad.filename}')

            # Create list of DFFFD flats to go in the DFFFD_flats stream
            elif 'FLAT' in ad.tags and ad.fiber_setup() in [DFFFD_setup]:
                flat_dfffd_list.append(ad)
                log.fullinfo(f'DFFFD Flat: {ad.filename}')

            # Warn if non-flats are in input - other fiber setup incorrect
            else:
                mislabeled.append(ad)
                log.warning('Not registered as Flat: %s', ad.filename)

        # TODO: raise an error when the flats streams splits DDDDF and DFFFD 
        # but one of them ends up empty.
        # Provide warnings if we do not have both types of flats
        if not flat_fdddf_list:
            log.warning('No FDDDF Flats in input list')
        if not flat_dfffd_list:
            log.warning('No DFFFD Flats in input list')

        self.streams['DFFFD_flats'] = flat_dfffd_list
        log.fullinfo(
            'Flat streams sotred in: p.streams["DFFFD_flats"] and p.streams["main"]'
        )
        return flat_fdddf_list

    def combineFlatStreams(self, adinputs=None, **params):
        """
        Combine two flat field streams into a single unified flat.

        This primitive takes two streams of pre-processed flat fields
        (typically 'FDDDF' flats in the main stream and 'DFFFD' flats in a
        secondary stream) and creates a combined flat field by taking the
        maximum value at each pixel position. This produces a flat field with
        all fibers illuminated (FFFFF flat configuration).

        Parameters
        ----------
        adinputs : list of AstroData
            Primary stream of flat objects (typically FDDDF flats)
        suffix : str
            Suffix to be added to output files
        stream_2 : str
            Name of the secondary stream to combine with the main stream
            (typically 'DFFFD_flats')

        Returns
        -------
        list of AstroData
            List containing a single AstroData object with the combined flat field
            data.
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))

        stream_2 = params['stream_2']

        if stream_2 not in self.streams.keys():
            # If stream_2 does not exist, return the provided input without modification
            log.fullinfo('Stream %s does not exist so nothing to transfer', stream_2)
            return adinputs

        # We expect both streams to have a length of 1 as we have a single DFFFD flat
        # and a single FDDDF flat after stacking. Log a warning otherwise.
        stream_length = len(self.streams[stream_2])
        adinputs_length = len(adinputs)

        if not adinputs_length == stream_length == 1:
            # Return the input without modification as we have unexpected stream lengths
            log.warning(
                'Unexpected stream lengths: %s and %s', adinputs_length, stream_length
            )
            return adinputs

        adout = deepcopy(adinputs[0])
        adinputs_2 = self.streams[stream_2]

        # Combine the data from the two streams by taking the max at each pixel
        adout[0].data = np.max([adinputs[0][0].data, adinputs_2[0][0].data], axis=0)
        # Do the same for the variance
        adout[0].variance = np.max(
            [adinputs[0][0].variance, adinputs_2[0][0].variance], axis=0
        )

        # For the rest of the extensions, we do not need to do this because we
        # will rerun id'ing on combined image frame
        return [adout]

    def splitBundle(self, adinputs=None, **params):
        """
        Split bundle into separate AstroData objects per extension.

        Split a bundle (Red and Blue multi-extension AstroData object) into
        separate AstroData objects, each containing one of the original extensions.

        This primitive creates a separate AstroData object for each extension
        in the input bundle. Each new object is a complete deep copy of the
        original, with all other extensions removed. This ensures all tables
        and associated data are preserved for the remaining extension. The
        original file name is retrieved from the 'ORIGNAME' header and assigned
        to the new object, and an 'ARCHNAME' header is added to reference the
        Gemini Archive file name. Additionally, certain header keywords in the
        EXPOSUREMETER table are renamed to avoid being lost by writing methods.

        Parameters
        ----------
        adinputs : list of AstroData objects
            List of multi-extension AstroData objects to be split.

        keep_suffix : bool
            If True, retains the suffix from the processed file name.

        Returns
        -------
        adouputs : list of AstroData objects
            List containing a separate AstroData object for each extension in each input
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        keep_suffix = params.get('keep_suffix', False)

        adoutputs = []

        for ad in adinputs:
            log.stdinfo(f"{ad.filename}: splitting bundle")
            for ext in ad.indices:
                # Create a deep copy of the original AstroData object
                arm_ad = deepcopy(ad)

                # Keep only the current extension and remove others
                indices_to_remove = [i for i in arm_ad.indices if i != ext]
                for i in indices_to_remove:
                    del arm_ad[i]

                # Revert to original file name
                arm_name = ad[ext].hdr.get('ORIGNAME')
                arm_ad.filename = arm_name

                # Overwrite ORIGNAME from PHU for consistency
                arm_ad.phu['ORIGNAME'] = arm_name

                # Now add a ARCHNAME for reference after ORIGNAME
                archive_card = fits.Card('ARCHNAME', ad.filename, 'Gemini archive name')
                arm_ad.phu.insert('ORIGNAME', archive_card, after=True)

                # Rename keywords of exposumeter table meta header
                # This is patch because the EXPOSUMETER.meta header looses
                # these particular keywords when writing to file
                if {'RAW', 'SCI'}.issubset(ad.tags):
                    arm_ad.EXPOSUREMETER.meta['header'].rename_keyword(
                        'TZERO2', 'ZP_PC'
                    )
                    arm_ad.EXPOSUREMETER.meta['header'].rename_keyword(
                        'TZERO3', 'ZP_FRD'
                    )

                # if processed file has sufix, add the same sufix to split files
                if 'PROCESSED' in ad.tags and keep_suffix:
                    suffix = ad.filename.split('_')[-1].replace('.fits', '')
                    arm_ad.update_filename(suffix='_' + suffix, strip=False)
                adoutputs.append(arm_ad)

        return adoutputs

    def fitDarkCoefficients(self, adinputs=None, **params):
        """
        Construct coefficients for log-linear fit of flux vs. exposure time.

        Parameters
        ----------
        adinputs : list of AstroData
            Input frames to be combined
        suffix : str
            Suffix to be added to output files

        Returns
        -------
        list of AstroData
            Combined output frame
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Validate inputs
        if len(adinputs) < 5:
            msg = f'Need at least 5 dark frames for fitting, got {len(adinputs)}'
            log.error(msg)
            raise ValueError(msg)

        # Extract data and metadata from AstroData objects
        framelist = []
        exposuretimes = []
        ndfilter = []
        filenames = []

        for ad in adinputs:
            # Get the data as float32
            frame = ad[0].data.astype(np.float32)
            framelist.append(frame)

            # Get exposure time
            exptime = ad.exposure_time()
            exposuretimes.append(float(exptime))

            # Get ND filter position
            nd_pos = float(ad.filter_orientation()['ND'])
            ndfilter.append(nd_pos)

            filenames.append(ad.filename)

        # Sort by exposure time (maintain correspondence with indices)
        sorted_indices = sorted(
            range(len(exposuretimes)), key=lambda i: exposuretimes[i]
        )

        framelist = [framelist[i] for i in sorted_indices]
        exposuretimes = [exposuretimes[i] for i in sorted_indices]
        ndfilter = [ndfilter[i] for i in sorted_indices]
        filenames = [filenames[i] for i in sorted_indices]

        log.fullinfo('Input files sorted by exposure time:')
        for filename, exptime in zip(filenames, exposuretimes):
            log.fullinfo('  %s: %ss', filename, exptime)

        # Create data cube
        log.fullinfo('Creating data cube from individual frames')
        cube = np.dstack(tuple(framelist))

        # Convert to log space for fitting
        exposuretimes = np.array(exposuretimes)
        logexptime = np.log10(exposuretimes)

        log.fullinfo('Fitting log-linear coefficients')
        # Initialize coefficient arrays
        z0 = np.empty(cube.shape[0:2], dtype=float)
        z1 = np.empty(cube.shape[0:2], dtype=float)

        # Fit polynomials
        for x in range(z0.shape[0]):
            for y in range(z0.shape[1]):
                z = np.polyfit(logexptime, cube[x, y, :], 1)
                z0[x, y] = z[0]
                z1[x, y] = z[1]

        log.fullinfo('Polynomial fitting completed')

        # Create output AstroData object based on first input
        ad_out = deepcopy(adinputs[0])

        # Main data array is emptied
        ad_out[0].data = np.zeros((1, 1))

        # Store coefficient arrays as extensions
        ad_out[0].COEFF_Z0 = z0
        ad_out[0].COEFF_Z1 = z1

        # Create table for exposure time information
        exptime_table = Table()
        exptime_table['logexptime'] = logexptime
        exptime_table['exptime'] = exposuretimes
        exptime_table['filename'] = filenames
        ad_out[0].LOGEXPTIME = exptime_table

        # Update metadata
        ad_out.phu.set(
            'NCOMBINE', len(adinputs), 'Number of darks used for coefficients'
        )

        # Timestamp and update filename
        gt.mark_history(ad_out, primname=self.myself(), keyword=timestamp_key)
        ad_out.update_filename(suffix=params['suffix'], strip=False)

        log.fullinfo(f'Dark coefficient arrays written to {ad_out.filename}')

        return [ad_out]


##############################################################################
# Below are the helper functions for the primitives in this module           #
##############################################################################


def _scaleCube(cube, scale_mode='first_frame'):
    """
    Scale data cube by normalizing each frame to a scale factor.

    The scaling can be done relative to the first frame or the
    mean of all frames.

    Parameters
    ----------
    cube : numpy.ndarray
        3D array of shape (height, width, n_frames) containing the frames to combine
    scale_mode : str
        Mode of scaling. Options are 'first_frame' or 'mean_frame'.
        - 'first_frame': Scale each frame to match the total flux of the first frame.
        - 'mean_frame': Scale each frame to match the mean total flux of all frames.

    Returns
    -------
    numpy.ndarray
        Scaled data cube with the same shape as the input
    """
    # Determine reference value based on scale_mode
    if scale_mode == 'first_frame':
        for i in range(cube.shape[-1]):
            a = np.sum(cube[:, :, i])
            cube[:, :, i] = cube[:, :, i] / a * np.sum(cube[:, :, 0])

    elif scale_mode == 'mean_frame':
        sums = []
        for i in range(cube.shape[-1]):
            a = np.sum(cube[:, :, i])
            sums = np.append(sums, a)

        total = np.sum(sums) / cube.shape[-1]
        for i in range(cube.shape[-1]):
            cube[:, :, i] = cube[:, :, i] / sums[i] * total

    else:
        msg = f'Unknown scale_mode: {scale_mode}.'
        raise ValueError(msg)

    return cube


def make_report_backgroundfit(ad, **kwargs):
    """
    Generate a PDF diagnostic report for straylight removal.

    Creates PDF with diagnostic plots showing each stage of the
    straylight removal process.

    Parameters
    ----------
    ad : AstroData
        The AstroData object being processed (used for filename).
    data_original : ndarray
        Original frame data before straylight removal.
    data_masked : ndarray
        Frame data with orders masked for background fitting.
    bkg1 : Background2D
        First iteration background model from photutils.
    bkg2 : Background2D
        Second iteration (correction) background model.
    result : ndarray
        Final straylight-subtracted data.
    max_val_masked : float
        Maximum value in the masked data used for scaling plots.

    Returns
    -------
    str
        Filename of the generated PDF report.
    """
    pdf_name = f"{ad.filename.replace('.fits', '_backgroundfit.pdf')}"

    # get array from kwargs
    data_original = kwargs.get('data_original')
    data_masked = kwargs.get('data_masked')
    bkg1 = kwargs.get('bkg1')
    bkg2 = kwargs.get('bkg2')
    result = kwargs.get('result')
    max_val_masked = kwargs.get('max_val_masked')

    # for legacy patch
    if hasattr(bkg1, 'background'):
        bkg1_background = bkg1.background
        bkg1_mesh_nmasked = bkg1.mesh_nmasked
        bkg2_background = bkg2.background
        box_size = bkg1.box_size[0]
    else:    
        bkg1_background = bkg1
        bkg1_mesh_nmasked = kwargs.get('bkg1_mesh_nmasked')
        bkg2_background = bkg2
        box_size = kwargs.get('box_size')
        pdf_name = f"{ad.filename.replace('.fits', '_backgroundfit_legacy.pdf')}"    
    

    with PdfPages(pdf_name) as pdf:
        figsize = (10, 8)
        # Figure 1: Raw frame, bias corrected
        fig1 = plt.figure(figsize=figsize)
        plt.title(
            f'Raw frame, bias corrected (max signal level: {int(max_val_masked):5d})'
        )
        plt.imshow(data_original, origin='lower', vmin=0, vmax=np.nanmax(data_masked))
        plt.colorbar()
        pdf.savefig(fig1)
        plt.close(fig1)

        # Figure 2: Raw frame, orders masked
        fig2 = plt.figure(figsize=figsize)
        plt.title('Raw frame, orders masked')
        plt.imshow(data_masked, origin='lower', vmin=0, vmax=np.nanmax(data_masked))
        plt.colorbar()
        pdf.savefig(fig2)
        plt.close(fig2)

        # Figure 3: Background model (first iteration)
        fig3 = plt.figure(figsize=figsize)
        plt.title('Background model')
        plt.imshow(bkg1_background, origin='lower', vmin=0, vmax=np.nanmax(data_masked))
        plt.colorbar()
        pdf.savefig(fig3)
        plt.close(fig3)

        # Figure 4: Background subtracted data (first iteration)
        fig4 = plt.figure(figsize=figsize)
        plt.title('Background subtracted data')
        plt.imshow(data_original - bkg1_background, origin='lower', vmin=-30, vmax=30)
        plt.colorbar()
        pdf.savefig(fig4)
        plt.close(fig4)

        # Figure 5: Number of masked pixels in background mesh
        fig5 = plt.figure(figsize=figsize)
        plt.title('Number of masked pixels in background mesh')
        #box_size = bkg1.box_size[0]
        plt.imshow(bkg1_mesh_nmasked, origin='lower', vmin=0, vmax=box_size**2)
        plt.colorbar()
        pdf.savefig(fig5)
        plt.close(fig5)

        # Figure 6: Background model, correction step
        fig6 = plt.figure(figsize=figsize)
        plt.title('Background model, correction step')
        plt.imshow(bkg2_background, origin='lower', vmin=-10, vmax=10)
        plt.colorbar()
        pdf.savefig(fig6)
        plt.close(fig6)

        # Figure 7: Final background subtracted data
        fig7 = plt.figure(figsize=figsize)
        plt.title('Final background subtracted data')
        plt.imshow(result, origin='lower', vmin=-30, vmax=30)
        plt.colorbar()
        pdf.savefig(fig7)
        plt.close(fig7)

    return pdf_name
