

from gempy.library import config
from geminidr.core import parameters_calibdb



class getProcessedDarkCoeffConfig(config.Config):
    pass


class storeCalibrationConfig(parameters_calibdb.storeCalibrationConfig):
    caltype = config.ChoiceField("Type of calibration", str,
                                 allowed={"processed_arc": "processed ARC",
                                          "processed_bpm": "processed BPM",
                                          "processed_dark": "processed DARK",
                                          "processed_dark_coeff": "processed DARK_COEFF",
                                          "processed_flat": "processed FLAT",
                                          "processed_sci": "processed SCI"},
                                 optional=False)


class storeProcessedDarkCoeffConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_darkCoefficients", optional=True)
