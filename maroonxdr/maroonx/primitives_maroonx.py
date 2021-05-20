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

    def subtractOverscan(self, adinputs=None, **params):
        """
        This primitive subtracts the overscan level from the image. The
        level for each row (currently the primitive requires that the overscan
        region be a vertical strip) is determined in one of the following
        ways, according to the *function* and *order* parameters:
        "poly":   a polynomial of degree *order* (1=linear, etc)
        "spline": using *order* equally-sized cubic spline pieces or, if
                  order=None or 0, a spline that provides a reduced chi^2=1
        "none":   no function is fit, and the value for each row is determined
                  by the overscan pixels in that row
        The fitting is done iteratively but, in the first instance, a running
        median of the rows is calculated and rows that deviate from this median
        are rejected (and used in place of the actual value if function="none")
        The GMOS-specific version of this primitive sets the "nbiascontam" and
        "order" parameters to their Gemini-IRAF defaults if they are None. It
        also removes the bottom 48 (ubinned) rows of the Hamamatsu CCDs from
        consideration in a polynomial fit. It then calls the generic version
        of the primitive.
        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        niterate: int
            number of rejection iterations
        high_reject: float
            number of standard deviations above which to reject high pixels
        low_reject: float
            number of standard deviations above which to reject low pixels
        overscan_section: str/None
            comma-separated list of IRAF-style overscan sections
        nbiascontam: int/None
            number of columns adjacent to the illuminated region to reject
        function: str
            function to fit ("polynomial" | "spline" | "none")
        order: int
            order of Chebyshev fit or spline/None
        """
        # To avoid crashing at the first line
        if not adinputs:
            return adinputs
        for ad in adinputs:
            dsec_list = ad.data_section()
            osec_list = ad.overscan_section()
            for ext, dsec, osec in zip(ad, dsec_list, osec_list):
                ext.hdr['BIASSEC'] = '[{}:{},{}:{}]'.format(osec.x1 + 1,
                                                            osec.x2, osec.y1 + 1, osec.y2)
                ext.hdr['DATASEC'] = '[{}:{},{}:{}]'.format(dsec.x1+1,
                                               dsec.x2, dsec.y1+1, dsec.y2)

        adinputs = super().subtractOverscan(adinputs, **params)
        return adinputs

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

