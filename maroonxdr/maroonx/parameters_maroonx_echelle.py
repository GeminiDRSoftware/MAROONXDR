# This parameter file contains the parameters related to the primitives
# defined in the primitive_maroonx_echelle.py file.

from gempy.library import config
from geminidr.core import parameters_stack
#from gemini.core import parameters_spect  # import core pkgs as needed.

class findStripesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")  # ?

