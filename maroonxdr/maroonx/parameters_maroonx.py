# This parameter file contains the parameters related to the primitives
# defined in the primitives_maroonx.py file, in alphabetical order

from gempy.library import config
from geminidr.core import parameters_ccd
from geminidr.core import parameters_stack, parameters_standardize

class addDQMXConfig(parameters_standardize.addDQConfig):
    suffix = config.Field("Filename suffix", str, "")
class checkArmConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class checkNDConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class combineFlatStreamsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")  # ?
    source = config.Field("Stream to transfer from", str, None)

class correctImageOrientationConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
class defineFlatStripesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
    slit_height = config.Field("Pixel illumination in cross-dispersion",
                               int, 10)
    extract = config.Field("Save extracted Stripes?", bool, False)

class findStripesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")

class identifyStripesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
    selected_fibers = config.Field("Fiber selection", str, None)

class overscanCorrectConfig(parameters_ccd.overscanCorrectConfig):
    def setDefaults(self):
        self.function = "none"

class removeStrayLightConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")
    box_size = config.Field("Photutils box size", int, 20)
    filter_size = config.Field("Photutilz filter size", int, 20)
class separateFlatStreamsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "")  # ?

class stackFramesMXCalConfig(parameters_stack.stackFramesConfig):
    suffix = config.Field("Filename suffix", str, "_stack")
    def setDefaults(self):
        self.scale = True
        self.zero = False
        self.separate_ext = True

class stackDarksConfig(parameters_stack.stackDarksConfig, stackFramesMXCalConfig):
    def setDefaults(self):
        self.reject_method = "sigclip"
        self.operation = "mean"
        self.hsigma = 2.  # dark stack is 2 sig, flat stack is 3 sig
        self.lsigma = 2.
        self.max_iters = 5
        self.scale = True

class stackFlatsConfig(parameters_stack.stackFlatsConfig,stackFramesMXCalConfig):
    suffix = config.Field("Filename suffix", str, "")  # ?
    def setDefaults(self):
        self.reject_method = "sigclip"
        self.operation = "mean"
        self.hsigma = 3.
        self.lsigma = 3.
        self.max_iters = 5
        self.scale = True

class validateDataConfig(parameters_standardize.validateDataConfig):
    suffix = config.Field("Filename suffix", str, "")
