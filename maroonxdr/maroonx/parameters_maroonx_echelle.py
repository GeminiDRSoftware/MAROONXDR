'''
This parameter file contains the parameters related to the
primitives defined in the primitive_maroonx_echelle.py file.
'''
from astrodata import AstroData

from gempy.library import config
from geminidr.core import parameters_ccd


class attachSyntheticDarkConfig(config.Config):
    '''
    This parameter set controls the attachDark primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")
    dark_coeff = config.ListField("Dark coefficient file", (str, AstroData), None,
                                optional=True, single=True)
    individual = config.Field("Unique dark for each frame", bool, False)

class attachDarkSubtractionConfig(config.Config):
    '''
    This parameter set controls the darkSubtraction primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")
    dark_type = config.ChoiceField("Dark to attach", str,
                                allowed={"synthetic": "Interpolated dark",
                                         "closest": "Closest in time"},
                                default="synthetic")

class createSyntheticDarkConfig(config.Config):
    '''
    This parameter set controls the createSyntheticDark primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")
    dark_coeff = config.ListField("Dark coefficient file", (str, AstroData), None,
                                optional=True, single=True)
    individual = config.Field("Unique dark for each frame", bool, False)

class extractStripesConfig(config.Config):
    '''
    This parameter set controls the extractStripes primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")
    skip_dark = config.ListField(
        'Input fibers for which the dark frame will NOT be subtracted.',
        int,
        None,
        optional=True,
        single=True,
    )
    remove_straylight = config.ListField(
        'Input fibers for which straylight will be removed.',
        int,
        None,
        optional=True,
        single=True,
    )
    slit_height = config.Field("Pixel illumination in cross-dispersion",int, 10)
    test_extraction = config.Field("Save in FITS-readable format for testing", bool, False)

class optimalExtractionConfig(config.Config):
    '''
    This parameter set controls the optimalExtraction primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "_reduced")
    opt_extraction = config.ListField("Fibers for optimal extraction", int, None, optional=True, single=False)
    full_output = config.Field("More outputs made", bool, False)
    penalty = config.Field("scaling penalty factor", float, None)
    s_clip = config.Field("sigma-clipping factor", float, None)
    back_var = config.Field("background variance", float, None)
    read_noise = config.Field("read noise", float, None)
    gain = config.Field("gain",float,None)
    def setDefaults(self):
        self.penalty = 1.0
        self.s_clip = 5.0
        self.back_var = 0.0
        self.read_noise = 1.14
        self.gain = 2.72

class boxExtractionConfig(config.Config):
    '''
    This parameter set controls the boxExtraction primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "")

# class darkSubtraction_oldConfig(config.Config):
#     '''
#     This parameter set controls the darkSubtraction primitive for MAROON-X.
#     '''
#     suffix = config.Field("Filename suffix", str, "")
#     individual = config.Field("individual or group caldb call", bool, False)

class overscanCorrectConfig(parameters_ccd.overscanCorrectConfig):
    '''
    This parameter set controls the overscanCorrect primitive for MAROON-X.
    '''
    def setDefaults(self):
        self.function = "none"
