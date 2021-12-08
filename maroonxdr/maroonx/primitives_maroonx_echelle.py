#
#                                                                       DRAGONS
#
#                                                 primitives_maroonx_echelle.py
# ------------------------------------------------------------------------------

from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.gemini import gemini_tools as gt
import numpy as np
from scipy import ndimage

from geminidr.core import Spect
from .primitives_maroonx import MAROONX
from . import parameters_maroonx_echelle

from recipe_system.utils.decorators import parameter_override
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
