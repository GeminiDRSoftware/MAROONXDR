"""
Container class for the extracted spectra of all fibers of an exposure.

``MXSpectrum`` reads the per-fiber extensions of a reduced AstroData
object and instantiates the matching spectrum class for each fiber
(``EchelleSpectrum``, ``EtalonSpectrum``, or ``FlatSpectrum``).
"""
from gempy.utils import logutils

from .etalonspectrum import EtalonSpectrum
from .flatspectrum import FlatSpectrum
from .echellespectrum import EchelleSpectrum


class MXSpectrum:
    """
    Extracted spectra of all fibers of a reduced MaroonX exposure.

    Reads the ``PEAKS`` table and the per-fiber extensions of the input
    AstroData object and instantiates one spectrum object per fiber:
    ``EtalonSpectrum`` for etalon fibers, ``FlatSpectrum`` for flat
    fibers, and ``EchelleSpectrum`` otherwise. Dark fibers and fibers
    without extracted data are skipped.

    Parameters
    ----------
    adinput : AstroData
        Reduced AstroData object with per-fiber box and optimal
        extraction, reduced orders, wavelength, and ``PEAKS``
        extensions.

    pm : PeakModeller
        Fit model for peaks. Inherited from the legacy pipeline;
        forwarded to the spectrum classes but currently unused.

    etalon_peaks_symmetric : bool
        If True, the etalon peaks were fit with equal left and right
        sigmas. Currently only logged; not forwarded to the spectrum
        classes. Default is False.

    wave_ext : str
        Name prefix of the wavelength solution extensions to load
        (``{wave_ext}_FIBER_{n}``). Default is 'WLS_STATIC'; the
        dynamic wavelength solution primitives pass 'WLS_DYNAMIC'.

    Attributes
    ----------
    spectra : dict
        Spectrum object per fiber number (1-5), or None for skipped
        fibers.

    echellogram : None
        Placeholder, currently unused.
    """
    def __init__(self, adinput, pm=None, etalon_peaks_symmetric=False, wave_ext='WLS_STATIC'):
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

