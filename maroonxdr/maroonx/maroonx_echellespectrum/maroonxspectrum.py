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
            the wavelength solution is shifted by this amount
        etalon_peaks_symmetric: bool
            if True, the etalon peaks are assumed to be symmetric around the central peak
        """
        self.logger = logutils.get_logger(__name__)
        logger = self.logger
        if etalon_peaks_symmetric:
            logger.utils("Using symmetric etalon peaks")

        # Check the fibers
        fibers = adinput.fiber_setup()

        poly_data = adinput[0].POLY
        peak_data = adinput[0].PEAKS
        # Convert peak_data to pandas dataframe
        peak_df = peak_data.to_pandas()
        
        self.spectra = {}
        self.echellogram = None
        '''
        Access data from input file.
        ROHAN:  This is where I stopped my work.
        I actually don't know if anything below this comment is useful.
        I think all we should need is the peaks per order
        And then we would apply a 30 knot spline to the peaks to get the wavelength solution
        '''

        for fiber_number, fiber in enumerate(fibers, start=1):

            reduced_orders = getattr(adinput[0], f'REDUCED_ORDERS_FIBER_{fiber_number}')
            box_data = getattr(adinput[0], f'BOX_REDUCED_FIBER_{fiber_number}')
            peak_data = peak_df.loc[peak_df['FIBER'] == fiber_number]
            
            if fiber == 'Etalon':
                # Create the EtalonSpectrum object
                self.spectra[fiber_number] = EtalonSpectrum(
                    box_data = box_data,
                    orders = reduced_orders,
                    peak_data = peak_data,
                    poly_data = poly_data,
                    pm = pm,
                    etalon_peaks_symmetric = etalon_peaks_symmetric)
            elif fiber == 'Flat':
                # Create the FlatSpectrum object
                self.spectra[fiber_number] = FlatSpectrum(box_data = box_data,
                                                          orders = reduced_orders,
                                                          peak_data = peak_data,
                                                          poly_data = poly_data,
                                                          pm = pm)

            elif fiber == 'Dark':
                self.spectra[fiber_number] = None

            else:
                # Treat as regular Echelle spectrum
                self.spectra[fiber_number] = EchelleSpectrum(box_data = box_data,
                                                             orders=reduced_orders,
                                                             peak_data = peak_data,
                                                             pm = pm)

