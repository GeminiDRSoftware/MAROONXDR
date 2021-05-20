# This parameter file contains the paramters related to the primitives
# define in the primitives_maroonx.py file

from gempy.library import config
from geminidr.core import parameters_ccd, parameters_nearIR, parameters_preprocess

class checkArmConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class checkNDConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class somePrimitiveConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_suffix")
    param1 = config.Field("Param1", str, "default")
    param2 = config.Field("do param2?", bool, False)

class someStuffConfig(config.Config):
    suffix = config.Field("Output suffix", str, "_somestuff")

class subtractOverscanConfig(parameters_ccd.subtractOverscanConfig):
    def setDefaults(self):
        self.function = "none"

    def validate(self):
        config.Config.validate(self)
        if self.function == "spline" and self.order == 0:
            raise ValueError("Must specify a positive spline order, or None")

class overscanCorrectConfig(subtractOverscanConfig, parameters_ccd.trimOverscanConfig):
    def setDefaults(self):
        self.suffix = "_overscanCorrected"

class stackDarks(parameters_nearIR.stackDarksConfig):
    def setDefaults(self):
        self.reject_method = "sigclip"
        self.hsigma = 2.
        self.lsigma = 2.