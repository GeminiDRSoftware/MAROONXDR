"""
This parameter file contains the parameters related to the primitives
defined in the primitives_maroonx_2D.py file, in alphabetical order
"""

from geminidr.core import parameters_ccd, parameters_stack, parameters_standardize
from gempy.library import config


class addDQMXConfig(parameters_standardize.addDQConfig):
    """
    This parameter set controls the addDQ primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')


class checkArmConfig(config.Config):
    """
    This parameter set controls the checkArm primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')


class checkMasterConfig(config.Config):
    """
    This parameter set controls the checkMaster primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')

class checkNDConfig(config.Config):
    """
    This parameter set controls the checkND primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')


class combineFlatStreamsConfig(config.Config):
    """
    This parameter set controls the combineFlatStreams primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')  # ?
    stream_2 = config.Field('Second stream to combine', str, None)


class correctImageOrientationConfig(config.Config):
    """
    This parameter set controls the correctImageOrientation primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')


class defineFlatStripesConfig(config.Config):
    """
    This parameter set controls the defineFlatStripes primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')
    slit_height = config.Field('Pixel illumination in cross-dispersion', int, 10)
    extract = config.Field('Save extracted Stripes?', bool, False)


class findStripesConfig(config.Config):
    """
    This parameter set controls the findStripes primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')


class identifyStripesConfig(config.Config):
    """
    This parameter set controls the identifyStripes primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')
    selected_fibers = config.Field('Fiber selection', str, None)


class overscanCorrectConfig(parameters_ccd.overscanCorrectConfig):
    """
    This parameter set controls the overscanCorrect primitive for MAROON-X.
    """

    def setDefaults(self):
        self.function = 'none'


class removeStrayLightConfig(config.Config):
    """
    This parameter set controls the removeStrayLight primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')
    snapshot = config.Field('save difference', bool, False)
    box_size = config.Field('Photutils box size', int, 20)
    filter_size = config.Field('Photutils filter size', int, 19)


class separateFlatStreamsConfig(config.Config):
    """
    This parameter set controls the separateFlatStreams primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')  # ?


class stackFramesMXCalConfig(parameters_stack.stackFramesConfig):
    """
    This parameter set controls the stackFrames primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '_stack')

    def setDefaults(self):
        self.scale = True
        self.zero = False
        self.separate_ext = True


class stackDarksConfig(parameters_stack.stackDarksConfig):
    """
    This parameter set controls the stackDarks primitive for MAROON-X.
    """

    # scale_mode can be 'first_frame' (darks) or 'mean_frame' (flats)
    scale_mode = config.Field('Scaling method for frames', str, 'first_frame')
    
    def setDefaults(self):
        self.reject_method = 'sigclip'
        self.hsigma = 2.0
        self.lsigma = 2.0
        self.max_iters = 5


class stackFlatsOldConfig(parameters_stack.stackFlatsConfig, stackFramesMXCalConfig):
    """
    This parameter set controls the stackFlats primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')

    def setDefaults(self):
        self.reject_method = 'sigclip'
        self.operation = 'mean'
        self.hsigma = 3.0
        self.lsigma = 3.0
        self.max_iters = 5
        self.scale = True

class stackFlatsConfig(parameters_stack.stackFlatsConfig, stackFramesMXCalConfig):
    """
    This parameter set controls the stackFlats primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')
    stream = config.Field('Stream name of flats to combine', str, 'main')
    # scale_mode can be 'first_frame' (darks) or 'mean_frame' (flats)
    scale_mode = config.Field('Scaling method for frames', str, 'mean_frame')

    def setDefaults(self):
        self.reject_method = 'sigclip'
        self.hsigma = 3.0
        self.lsigma = 3.0
        self.max_iters = 5

class validateDataConfig(parameters_standardize.validateDataConfig):
    """
    This parameter set controls the validateData primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')

    def setDefaults(self):
        self.require_wcs = False


class splitBundleConfig(config.Config):
    """
    This parameter set controls the splitBundle primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')


class fitDarkCoefficientsConfig(config.Config):
    """
    This parameter set controls the fitDarkCoefficients primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')