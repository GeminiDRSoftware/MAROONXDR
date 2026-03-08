__all__ = ['AstroDataMAROONX']

from astrodata import factory
from gemini_instruments.gemini import addInstrumentFilterWavelengths
from .adclass import AstroDataMAROONX
from .lookup import filter_wavelengths

factory.addClass(AstroDataMAROONX)

addInstrumentFilterWavelengths('MAROONX', filter_wavelengths)

# Register MAROONX calibration association rules with FitsStorage
from fits_storage.cal.calibration import inst_class
from .calibration_maroonx import CalibrationMAROONX
inst_class["MAROONX"] = CalibrationMAROONX

# Register MaroonX-specific caltypes with the local calibration manager
from recipe_system.cal_service.localmanager import args_for_cals
args_for_cals['processed_wavecal'] = ('wavecal', {'processed': True})
args_for_cals['processed_dark_coeff'] = ('dark_coeff', {'processed': True})

