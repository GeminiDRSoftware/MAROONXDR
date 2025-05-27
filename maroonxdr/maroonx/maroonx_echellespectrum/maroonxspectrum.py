from gempy.utils import logutils
import numpy as np
from .etalonspectrum import EtalonSpectrum
from .flatspectrum import FlatSpectrum
from .echellespectrum import EchelleSpectrum

def save_wavelength_solution(wavelength_solution, filename, ext_name = "wavelength_solution", overwrite=False):
    """
    Saves the wavelength solution to a fits file in a new extension.

    Args:
        wavelength_solution : Wavelength solution to save.
        filename (str): Name of the file to save the wavelength solution to.
        ext_name (str): Name of the extension to save the wavelength solution to.
        overwrite (bool): Overwrite extension if it already exists.

    Returns:
        None
    """
    pass

def load_wavelength_solution(adinput, ext_name = "wavelength_solution"):
    """
    Loads the wavelength solution from a fits file.

    Args:
        filename (str): Name of the file to load the wavelength solution from.
        ext_name (str): Name of the extension to load the wavelength solution from.

    Returns:
        wavelength_solution : Wavelength solution loaded from the file.
    """
    pass

class MXSpectrum:
    '''
    This class is used to read in a MaroonX spectrum and apply the wavelength solution.
    '''
    def __init__(self, adinput, pm=None, etalon_peaks_symmetric=False):
        """
        Initializes the MXSpectrum object.

        Parameters
        ----------
        adinput: AstroData object
            the AstroData object to process
        pm: float
            Fit model for peaks. If given, it will be used when fitting lines. 
            Otherwise a Gaussian model will be used.
        etalon_peaks_symmetric: bool
            if True, the etalon peaks are assumed to be symmetric around the central peak
        """
        self.logger = logutils.get_logger(__name__)
        logger = self.logger
        if etalon_peaks_symmetric:
            logger.utils("Using symmetric etalon peaks")

        # Check the fibers
        fibers = adinput.fiber_setup()

        # poly_data = adinput[0].POLY
        peak_data = adinput[0].PEAKS.to_pandas()
        peak_data = peak_data.sort_values(by=['FIBER', 'ORDER', 'CENTER'])
        
        self.spectra = {}
        self.echellogram = None
        '''
        Access data from input file.
        ROHAN:  This is where I stopped my work.
        I actually don't know if anything below this comment is useful.
        I think all we should need is the peaks per order
        And then we would apply a 30 knot spline to the peaks to get the wavelength solution
        '''

        # Define the spectra classes based on fiber type
        spectra_classes = {
            'Echelle': EchelleSpectrum,
            'Etalon': EtalonSpectrum,
            'Flat lamp': FlatSpectrum,
        }

        for fiber_number, fiber in enumerate(fibers, start=1):

            if fiber == 'Dark':
                # Skip Dark fiber
                self.spectra[fiber_number] = None
                continue

            reduced_orders = getattr(adinput[0], f'REDUCED_ORDERS_FIBER_{fiber_number}')
            box_data = getattr(adinput[0], f'BOX_REDUCED_FIBER_{fiber_number}')
            peaks = peak_data.loc[peak_data['FIBER'] == fiber_number]
            wls_static_data = getattr(adinput[0], f'WLS_STATIC_FIBER_{fiber_number}')

            # If the fiber is not in the spectra_classes, default to EchelleSpectrum
            spectra_cls = spectra_classes.get(fiber, EchelleSpectrum)
            self.spectra[fiber_number] = spectra_cls(
                box_data=box_data,
                orders=reduced_orders,
                peak_data=peaks,
                wavelength_data=wls_static_data,
                fiber=fiber_number,
                pm=pm,
            )

