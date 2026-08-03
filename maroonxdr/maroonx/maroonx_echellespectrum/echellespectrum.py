"""
Base class for extracted 1-D echelle spectra.

``EchelleSpectrum`` holds the extracted data of a single fiber, one row
per echelle order, and applies wavelength solutions to it. All other
spectrum types in this package inherit from it, and it is used directly
for target fibers.
"""

import pandas as pd
import numpy as np

from gempy.gemini import gemini_tools as gt

from astropy.utils.decorators import deprecated

_DEPRECATION_MSG = (
    "This method is inherited from the legacy MaroonX pipeline and "
    "is not used by the DRAGONS reduction; it may be removed in a "
    "future release."
)

class EchelleSpectrum:
    r"""
    Extracted 1-D echelle spectrum of a single fiber.

    Holds the box and optimal extracted intensity data, the associated
    errors, and the wavelengths for each order. The per-order data is
    stored in a pandas DataFrame indexed by the sorted physical orders,
    exposed through the ``data`` property.

    Parameters
    ----------
    orders : array-like of int
        Physical echelle orders.

    box_data : ndarray
        Box extracted intensity data vector for each order. Shape
        (N_orders, N_pixels), where N_pixels is the number of pixels in
        the dispersion direction.

    box_error : ndarray
        Error vector for each order. For a CCD gain of 1, this is the
        square root of the intensity.

    opt_data : ndarray
        Optimal extracted intensity data vector for each order. Shape
        (N_orders, N_pixels).

    opt_error : ndarray
        Optimal extracted variance vector for each order. Shape
        (N_orders, N_pixels).

    wavelength_data : ndarray
        Wavelength data vector for each order (nm). Shape
        (N_orders, N_pixels).

    fiber : int
        Fiber number. Default is 1.

    pm : PeakModeller
        Fit model for peaks. Inherited from the legacy pipeline, where it
        was used when fitting ThAr lines; stored but currently unused.

    filename : str
        Filename of the raw file. Inherited from the legacy pipeline;
        accepted but not stored.

    \*\*kwargs
        Additional keyword arguments, silently ignored.

    Attributes
    ----------
    norm_orders : ndarray
        Orders normalized to the range [-1, 1].

    etalon_parameters : lmfit.parameter.Parameters or None
        Etalon parameters. None until set by ``EtalonSpectrum`` or by
        the wavelength solution primitives.

    model : None
        Placeholder, currently unused.
    """
    def __init__(self, orders, box_data=None, box_error=None
                 , opt_data=None, opt_error=None, wavelength_data=None, fiber=1, pm=None,
                 filename=None, **kwargs):
        # Convert orders to integers if needed
        orders = np.asarray(orders, dtype=int)
        
        self.box_data = box_data
        self.box_error = box_error
        self.opt_data = opt_data
        self.opt_error = opt_error
        self.wavelength_data = wavelength_data

        _tolist = lambda x: None if x is None else x.tolist()
        self._data = pd.DataFrame({
            'box_data': _tolist(box_data),
            'box_error': _tolist(box_error),
            'opt_data': _tolist(opt_data),
            'opt_error': _tolist(opt_error),
            'wavelength': _tolist(wavelength_data),
        }, index=orders)
        self._data.sort_index(axis=0, inplace=True)

        # physical orders - sorted like the index
        self._orders = np.sort(orders)
        min_order = np.min(orders)
        max_order = np.max(orders)
        # normalized orders
        self.norm_orders = self.normalize_orders(orders, min_order, max_order)

        self.etalon_parameters = None
        self.model = None

        self.fiber = fiber
        self.pm = pm

        self.log = gt.logutils.get_logger(__name__)

    @property
    def data(self):
        """Return the per-order data as a DataFrame indexed by physical order."""
        return self._data

    @property
    def orders(self):
        """Return the sorted physical echelle orders."""
        return self._orders

    def normalize_orders(self, order, min_order, max_order):
        """
        Convert physical orders to normalized orders in the range [-1, 1].

        Parameters
        ----------
        order : int or ndarray
            Physical order(s).

        min_order : int
            Minimum physical order. Mapped to -1.

        max_order : int
            Maximum physical order. Mapped to 1.

        Returns
        -------
        float or ndarray
            Normalized order(s).
        """
        norm_order = (order - min_order)/(max_order - min_order) * 2. - 1.
        return norm_order

    def normalize_pixel(self, pixel):
        """
        Convert physical pixels to normalized pixels in the range [-1, 1].

        The normalization is based on the number of pixels in the
        dispersion direction of ``box_data``.

        Parameters
        ----------
        pixel : int or ndarray
            Physical pixel(s).

        Returns
        -------
        float or ndarray
            Normalized pixel(s).
        """
        # norm_pixel = (pixel - min_pixel)/(max_pixel - min_pixel) * 2. - 1.
        norm_pixel = pixel / self.box_data.shape[1] * 2. - 1.
        return norm_pixel

    @deprecated(since="DRAGONS-integration", message=_DEPRECATION_MSG)
    def data_flattened(self, box_data=False):
        """
        Return the data flattened across orders.

        Parameters
        ----------
        box_data : bool
            If False, data is optimal extraction, if True, data is box
            extraction.

        Returns
        -------
        tuple of ndarray
            ``(intensity, wavelength)``, each flattened across orders.
        """
        data_selection = 'box_data' if box_data else 'opt_data'
        return np.hstack(self.data[data_selection].values), np.hstack(self.data['wavelength'].values)

    @deprecated(since="DRAGONS-integration", message=_DEPRECATION_MSG)
    def min_wavelength(self):
        """Return the minimum wavelength (nm)."""
        return np.min(self.data['wavelength'])

    @deprecated(since="DRAGONS-integration", message=_DEPRECATION_MSG)
    def max_wavelength(self):
        """Return the maximum wavelength (nm)."""
        return np.max(self.data['wavelength'])

    @deprecated(since="DRAGONS-integration", message=_DEPRECATION_MSG)
    def find_orders_containing_wavelength(self, wavelength):
        """
        Find the orders that contain the specified wavelength.

        TODO: Optimize this function.

        Parameters
        ----------
        wavelength : float
            Wavelength to find orders for (nm).

        Returns
        -------
        list of ndarray
            Wavelength arrays of the orders that contain the specified
            wavelength.
        """
        found = []
        for i in self.data['wavelength']:
            if min(i) < wavelength < max(i):
                found.append(i)
        return found

    @deprecated(since="DRAGONS-integration", message=_DEPRECATION_MSG)
    def blaze_correct(self, flat_spectrum, box_data = False):
        """
        Add deblazed values to the data for either box or optimal extraction.

        Parameters
        ----------
        flat_spectrum : MXSpectrum
            Flat spectrum containing a ``FlatSpectrum`` for this fiber.

        box_data : bool
            If False, data is optimal extraction, if True, data is box
            extraction.
        """
        log = self.log
        data_selection = 'box_data' if box_data else 'opt_data'
        try:
            if len(flat_spectrum.spectra[self.fiber].blaze) == 0:
                log.info(f'Creating blaze functions for fiber {self.fiber}')
                flat_spectrum.spectra[self.fiber].fit_blaze()
            deblazed_int = []
            for o in self.orders:
                deblazed_int.append(self.data.loc[o][data_selection] / flat_spectrum.spectra[self.fiber].blaze[o])
            self._data["blaze_corrected_"+data_selection] = deblazed_int
        except:
            log.error('Blaze correction failed.  Flat spectrum may not have been loaded.')

    def apply_wavelength_solution(self, wavelength_solution):
        """
        Apply the specified wavelength solution to the spectrum.

        The wavelength column of the per-order data is replaced in place
        with the solution evaluated at every pixel of each order.

        Parameters
        ----------
        wavelength_solution : WavelengthSolution
            Wavelength solution to apply.
        """
        log = self.log
        log.info('Applying wavelength solution')

        for o in self.orders:
            x = np.arange(len(self.data.loc[o]['box_data']))
            self._data.loc[o]["wavelength"] = wavelength_solution.get_wavelength(x, o)
