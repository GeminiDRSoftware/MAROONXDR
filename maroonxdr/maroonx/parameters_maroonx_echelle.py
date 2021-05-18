# This parameter file contains the parameters related to the primitives
# defined in the primitive_maroonx_echelle.py file.

from gempy.library import config
#from gemini.core import parameters_spect  # import core pkgs as needed.

class myNewPrimitive(config.Config):
    suffix = config.Field("Filename suffix", str, "_suffix")
    param1 = config.Field("Param1", str, "default")
    param2 = config.Field("do param2?", bool, False)

class correctImageOrientation(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class findStripes(config.Config):
    suffix = config.Field("Filename suffix", str, "")  # ?