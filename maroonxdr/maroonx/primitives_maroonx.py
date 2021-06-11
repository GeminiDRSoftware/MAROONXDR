#
#                                                                       DRAGONS
#
#                                                         primitives_maroonx.py
# ------------------------------------------------------------------------------

from gempy.gemini import gemini_tools as gt
import numpy as np
import copy
from geminidr.gemini.primitives_gemini import Gemini
from geminidr.core import CCD, NearIR, primitives_preprocess

from . import parameters_maroonx

from .lookups import timestamp_keywords as maroonx_stamps

from recipe_system.utils.decorators import parameter_override
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

    def correctImageOrientation(self, adinputs=None, debug_level=0):
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

        for ad in adinputs:

            if debug_level > 0:
                plt.figure()
                plt.title('Original Image')
                plt.imshow(ad.data[0], origin='lower')

            if ad.image_orientation()['vertical orientation flip']:  # flip up-down
                ad.data[0] = np.flipud(ad.data[0])
            if ad.image_orientation()['horizontal orientation flip']:  # flip left-right (dispersion direction)
                ad.data[0] = np.fliplr(ad.data[0])

            if debug_level > 0:
                plt.figure()
                plt.title('Orientation Corrected image')
                plt.imshow(ad.data[0], origin='lower')
                plt.show()

        # return imgs
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)

        return adinputs

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
        for ad in adinputs:
            if check_val != ad.filter_orientation()['ND']:
                log.warning("Not all frames have the same simcal ND filter setting, restricting set to first seen")
            else:
                ad.update_filename(suffix=params['suffix'], strip=True)
                adoutputs.append(ad)
        if len(adoutputs) == 1:
            raise IOError("Less than two frames found with first frame simcal ND filter setting")
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
                log.fullinfo("FDDDF Flat: {}, {}".format(ad.data_label(), ad.filename))
            elif "FLAT" in tags and ad.fiber_setup() == ['Dark', 'Flat', 'Flat', 'Flat', 'Dark']:
                flat_DFFFD_list.append(ad)
                log.fullinfo("DFFFD Flat: {}, {}".format(ad.data_label(), ad.filename))
            else:
                mislabeled.append(ad)
                log.warning("Not Flat: {} {}".format(ad.data_label(),
                                                     ad.filename))
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

        adout = copy.deepcopy(adinputs[0])
        # add header updates (all fibers set to 'Flat', etc)
        # adout[0] *= 0  # this is probably not how to do this    356
        # adout[0] += np.max([adinputs[0].data[0], self.streams[source][0].data[0]], axis=0)
        adout[0].data = np.max([adinputs[0].data[0], self.streams[source][0].data[0]], axis=0)

        fromad2 = np.where(adout[0].data == self.streams[source][0].data[0])

        fromad = np.where(adout[0].data == adinputs[0].data[0])

        listcoo2 = list(zip(fromad2[0], fromad2[1]))

        listcoo = list(zip(fromad[0], fromad[1]))

        for coo in listcoo2:
            adout[0].variance[coo[0], coo[1]] = self.streams[source][0].variance[coo[0], coo[1]]

        for coo in listcoo:
            adout[0].variance[coo[0], coo[1]] = adinputs[0].variance[coo[0], coo[1]]

        return adout


    @staticmethod
    def _has_valid_extensions(ad):
        """ Check that the AD has a valid number of extensions. """

        # this needs to be updated at appropriate. 
        return len(ad) in [1]

