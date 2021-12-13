# This parameter file contains the paramters related to the primitives
# define in the primitives_maroonx.py file

from gempy.library import config
from geminidr.core import parameters_ccd, parameters_nearIR, parameters_stack

class checkArmConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class checkNDConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class correctImageOrientationConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class subtractOverscanConfig(parameters_ccd.subtractOverscanConfig):
    def setDefaults(self):
        self.function = "none"

class StackFramesConfig(parameters_stack.stackFramesConfig):
    suffix = config.Field("Filename suffix", str, "")  # ?
    def setDefaults(self):
        self.separate_ext = False

class stackDarksConfig(parameters_nearIR.stackDarksConfig,StackFramesConfig):
    def setDefaults(self):
        self.reject_method = "sigclip"
        self.operation = "median"
        self.hsigma = 2.
        self.lsigma = 2.
        self.max_iters = 5

class separateFlatStreamsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")  # ?

class stackFlatsConfig(parameters_stack.stackFlatsConfig,StackFramesConfig):
    suffix = config.Field("Filename suffix", str, "")  # ?
    def setDefaults(self):
        self.reject_method = "sigclip"
        self.operation = "median"
        self.hsigma = 3.
        self.lsigma = 3.
        self.max_iters = 5

class combineFlatStreamsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")  # ?
    source = config.Field("Stream to transfer from", str, None)

class findStripesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class identifyStripesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
    selected_fibers = config.Field("Fiber selection", str, None)

class defineFlatStripesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
    slit_height = config.Field("Pixel illumination in cross-dispersion", int, 10)
    extract = config.Field("Save extracted Stripes?", bool, False)

class removeStrayLightConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
    box_size = config.Field("Photutils box size", int, 20)
    filter_size = config.Field("Photutilz filter size", int, 20)