from gempy.gemini import gemini_tools as gt

from geminidr.core import CalibDB
from . import parameters_calibdb_maroonx

from recipe_system.utils.decorators import parameter_override

# ------------------------------------------------------------------------------
@parameter_override
class CalibDBMaroonX(CalibDB):
    """
    This is the class containing all of the calibration bookkeeping primitives
    for the CalibDBMaroonX level of the type hierarchy tree. It inherits all
    the primitives from the level above
    """
    tagset = set()  # Not allowed to be a selected as a primitivesClass

    def __init__(self, adinputs, **kwargs):
        super(CalibDBMaroonX, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'maroonxdr.maroonx.lookups'
        self._param_update(parameters_calibdb_maroonx)


    def getProcessedWavecal(self, adinputs=None, **params):
        raise NotImplementedError("This primitives is not implemented")

    # =========================================================
    # STORE PRIMITIVES
    # =========================================================

    def storeProcessedWavecal(self, adinputs=None, suffix=None):
        caltype = 'processed_wavecal'
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))
        adinputs = self._markAsCalibration(adinputs, suffix=suffix,
                                    primname=self.myself(), keyword="PRWAVECAL")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs