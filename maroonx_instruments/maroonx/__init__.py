__all__ = ['AstroDataMAROONX']

from astrodata import factory
from gemini_instruments.gemini import addInstrumentFilterWavelengths
from .adclass import AstroDataMAROONX
from .lookup import filter_wavelengths

factory.addClass(AstroDataMAROONX)

addInstrumentFilterWavelengths('MAROONX', filter_wavelengths)

