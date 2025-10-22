from gempy.gemini import gemini_tools as gt

from geminidr.core import CalibDB
from . import parameters_calibdb_maroonx


from recipe_system.utils.decorators import parameter_override

# Extend REQUIRED_TAG_DICT with MaroonX-specific calibration types
from recipe_system.cal_service.caldb import REQUIRED_TAG_DICT
try:
    REQUIRED_TAG_DICT['processed_wavecal'] = ['PROCESSED', 'WAVECAL']
except (TypeError, AttributeError):
    # Handle mocked environment (e.g., during documentation builds)
    pass

# ------------------------------------------------------------------------------
@parameter_override
class CalibDBMAROONX(CalibDB):
    """
    This is the class containing all of the calibration bookkeeping primitives
    for the CalibDBMAROONX level of the type hierarchy tree. It inherits all
    the primitives from the level above
    """
    tagset = set()  # Not allowed to be a selected as a primitivesClass

    def _initialize(self, adinputs, **kwargs):
        self.inst_lookups = "maroonxdr.maroonx.lookups"
        super()._initialize(adinputs, **kwargs)
        self._param_update(parameters_calibdb_maroonx)


    def getProcessedWavecal(self, adinputs=None, **params):
        procmode = 'sq' if self.mode == 'sq' else None
        cals = self.caldb.get_processed_wavecal(adinputs, procmode=procmode)
        self._assert_calibrations(adinputs, cals)
        return adinputs

    # =========================== STORE PRIMITIVES =================================
    def storeProcessedWavecal(self, adinputs=None, suffix=None, force=False):
        caltype = 'processed_wavecal'
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))

        adinputs = self._markAsCalibration(adinputs, suffix=suffix, update_datalab=True,
                                    primname=self.myself(), keyword="PRWAVECAL")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs