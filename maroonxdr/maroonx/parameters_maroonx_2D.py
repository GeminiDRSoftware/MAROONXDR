'''
This parameter file contains the parameters related to the primitives
defined in the primitives_maroonx_2D.py file, in alphabetical order
'''
from gempy.library import config
from geminidr.core import parameters_ccd
from geminidr.core import parameters_stack, parameters_standardize

class addDQMXConfig(parameters_standardize.addDQConfig):
    '''
    This parameter set controls the addDQ primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")

class checkArmConfig(config.Config):
    '''
    This parameter set controls the checkArm primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")

class checkNDConfig(config.Config):
    '''
    This parameter set controls the checkND primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")

class combineFlatStreamsConfig(config.Config):
    '''
    This parameter set controls the combineFlatStreams primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")  # ?
    source = config.Field("Stream to transfer from", str, None)

class correctImageOrientationConfig(config.Config):
    '''
    This parameter set controls the correctImageOrientation primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")

class defineFlatStripesConfig(config.Config):
    '''
    This parameter set controls the defineFlatStripes primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")
    slit_height = config.Field("Pixel illumination in cross-dispersion",
                               int, 10)
    extract = config.Field("Save extracted Stripes?", bool, False)

class findStripesConfig(config.Config):
    '''
    This parameter set controls the findStripes primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")

class identifyStripesConfig(config.Config):
    '''
    This parameter set controls the identifyStripes primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")
    selected_fibers = config.Field("Fiber selection", str, None)

class overscanCorrectConfig(parameters_ccd.overscanCorrectConfig):
    '''
    This parameter set controls the overscanCorrect primitive for MAROON-X.
    '''
    def setDefaults(self):
        self.function = "none"

class removeStrayLightConfig(config.Config):
    '''
    This parameter set controls the removeStrayLight primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")
    snapshot = config.Field("save difference", bool, False)
    box_size = config.Field("Photutils box size", int, 20)
    filter_size = config.Field("Photutilz filter size", int, 20)

class separateFlatStreamsConfig(config.Config):
    '''
    This parameter set controls the separateFlatStreams primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")  # ?

class stackFramesMXCalConfig(parameters_stack.stackFramesConfig):
    '''
    This parameter set controls the stackFrames primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "_stack")
    def setDefaults(self):
        self.scale = True
        self.zero = False
        self.separate_ext = True

class stackDarksConfig(parameters_stack.stackDarksConfig, stackFramesMXCalConfig):
    '''
    This parameter set controls the stackDarks primitive for MAROON-X.
    '''
    def setDefaults(self):
        self.reject_method = "sigclip"
        self.operation = "mean"
        self.hsigma = 3.
        self.lsigma = 3.
        self.max_iters = 5

class stackFlatsConfig(parameters_stack.stackFlatsConfig,stackFramesMXCalConfig):
    '''
    This parameter set controls the stackFlats primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")
    def setDefaults(self):
        self.reject_method = "sigclip"
        self.operation = "mean"
        self.hsigma = 3.
        self.lsigma = 3.
        self.max_iters = 5
        self.scale = True

class validateDataConfig(parameters_standardize.validateDataConfig):
    '''
    This parameter set controls the validateData primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")
