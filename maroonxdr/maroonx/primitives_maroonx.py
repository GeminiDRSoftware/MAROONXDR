#
#                                                                       DRAGONS
#
#                                                         primitives_maroonx.py
# ------------------------------------------------------------------------------

from gempy.gemini import gemini_tools as gt

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

    def someStuff(self, adinputs=None, **params):
        """
        Write message to screen.  Test primitive.

        Parameters
        ----------
        adinputs
        params

        Returns
        -------

        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        for ad in adinputs:
            log.status('I see '+ad.filename)

            gt.mark_history(ad, primname=self.myself(), keyword="TEST")
            ad.update_filename(suffix=params['suffix'], strip=True)

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


    def separateFlats(self, adinputs=None, **params):
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


    @staticmethod
    def _has_valid_extensions(ad):
        """ Check that the AD has a valid number of extensions. """

        # this needs to be updated at appropriate. 
        return len(ad) in [1]

