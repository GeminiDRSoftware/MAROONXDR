# This parameter file contains the parameters related to the primitives
# defined in the primitive_maroonx_echelle.py file.

from gempy.library import config
from geminidr.core import parameters_stack
#from gemini.core import parameters_spect  # import core pkgs as needed.

class extractStripesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
    slit_height = config.Field("Pixel illumination in cross-dispersion",int, 10)

class optimalExtractionConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
    full_output = config.Field("More outputs made", bool, False)
    penalty = config.Field("scaling penalty factor", float, None)
    s_clip = config.Field("sigma-clipping factor", float, None)