from gempy.utils import logutils

from .etalonspectrum import EtalonSpectrum
from .flatspectrum import FlatSpectrum
from .echellespectrum import EchelleSpectrum


class MXSpectrum:
    '''
    This class is used to read in a MaroonX spectrum and apply the wavelength solution.
    '''
    def __init__(self, adinput, pm=None, etalon_peaks_symmetric=False, wave_ext='WLS_STATIC'):
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
        peak_data['ORDER'] = peak_data['ORDER'].map(int)
        peak_data = peak_data.sort_values(by=['FIBER', 'ORDER', 'CENTER'])
        peak_data = peak_data.set_index(['FIBER', 'ORDER', 'CENTER'], drop=False)

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
            'Target': EchelleSpectrum,
            'Etalon': EtalonSpectrum,
            'Flat lamp': FlatSpectrum,
        }

        for fiber_number, fiber in enumerate(fibers, start=1):

            if fiber == 'Dark':
                # Skip Dark fiber
                self.spectra[fiber_number] = None
                continue

            reduced_orders = getattr(adinput[0], f'REDUCED_ORDERS_FIBER_{fiber_number}', None)
            box_data = getattr(adinput[0], f'BOX_REDUCED_FIBER_{fiber_number}', None)
            box_var = getattr(adinput[0], f'BOX_REDUCED_VAR_{fiber_number}', None)
            opt_data = getattr(adinput[0], f'OPTIMAL_REDUCED_FIBER_{fiber_number}', None)
            opt_var = getattr(adinput[0], f'OPTIMAL_REDUCED_VAR_{fiber_number}', None)
            wls_data = getattr(adinput[0], f'{wave_ext}_FIBER_{fiber_number}', None)

            if reduced_orders.size == 1:
                logger.warning(f"Missing data for fiber {fiber_number}. Skipping.")
                self.spectra[fiber_number] = None
                continue

            if fiber != 'Target':
                peaks = peak_data.loc[fiber_number].copy()
            else:
                peaks = None

            # If the fiber is not in the spectra_classes, default to EchelleSpectrum
            spectra_cls = spectra_classes.get(fiber, EchelleSpectrum)
            self.spectra[fiber_number] = spectra_cls(
                peak_data=peaks,
                box_data=box_data,
                box_err=box_var,
                opt_data=opt_data,
                opt_err=opt_var,
                orders=reduced_orders,
                wavelength_data=wls_data,
                fiber=fiber_number,
                pm=pm,
            )

