# This parameter file contains the parameters related to the primitives
# defined in the primitive_maroonx_echelle.py file.

from gempy.library import config
from geminidr.core import parameters_stack, parameters_standardize
from geminidr.core import parameters_ccd, parameters_calibdb  # import core pkgs as needed.


class extractStripesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
    slit_height = config.Field("Pixel illumination in cross-dispersion",int, 10)

class optimalExtractionConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_reduced")
    full_output = config.Field("More outputs made", bool, False)
    penalty = config.Field("scaling penalty factor", float, None)
    s_clip = config.Field("sigma-clipping factor", float, None)
    back_var = config.Field("background variance", float, None)
    read_noise = config.Field("read noise", float, None)
    gain = config.Field("gain",float,None)
    def setDefaults(self):
        self.penalty = 1.0
        self.s_clip = 5.0
        self.back_var = 200.0
        self.read_noise = 1.14
        self.gain = 2.72

class overscanCorrectConfig(parameters_ccd.overscanCorrectConfig):
    def setDefaults(self):
        self.function = "none"

# class storeProcessedScience(parameters_calibdb.storeProcessedScienceConfig):
#     def setDefaults(self):
#         self.update_datalab
