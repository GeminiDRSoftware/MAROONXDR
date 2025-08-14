
from gempy.library import config
from geminidr.core import parameters_calibdb



class getProcessedWavecalConfig(config.Config):
    pass



class storeCalibrationConfig(parameters_calibdb.storeCalibrationConfig):
    caltype = config.ChoiceField("Type of calibration", str,
                                 allowed={"processed_wavecal": "processed WAVECAL",
                                          "processed_bias": "procsessed BIAS",
                                          "processed_bpm": "processed BPM",
                                          "processed_dark": "processed DARK",
                                          "processed_flat": "processed FLAT"},
                                 optional=False)


class storeProcessedWavecalConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_wavecal", optional=True)

