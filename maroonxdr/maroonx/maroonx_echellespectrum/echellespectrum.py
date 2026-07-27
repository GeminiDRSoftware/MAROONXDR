'''
This class contains the echelle spectrum, from which all other spectrum types inherit.
'''

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
    '''
    The echelle spectrum class contains all information about extracted 1-D echelle spectra.
    Each object contains several orders, and corresponding wavelength and intensity data.
    '''
    def __init__(self, orders, box_data=None, box_error=None
                 , opt_data=None, opt_error=None, wavelength_data=None, fiber=1, pm=None,
                 filename=None, **kwargs):
        """
        Loads an echelle spectrum.

        Parameters
        ----------
        orders : list
            List of orders.
        box_data : numpy array
            box extracted intensity data vector for each order. Shape (N_orders, N_pixels), where N_pixels is the
            number of pixels in the dispersion direction.
        box_error : numpy array
            error vector for intensity order.  For a CCD gain of 1, this is the square root of the intensity.
        opt_data : numpy array
            optimal extracted intensity data vector for each order. Shape (N_orders, N_pixels), where N_pixels is the
            number of pixels in the dispersion direction.
        opt_error : numpy array
            optimal extracted variance vector for each order. Shape (N_orders, N_pixels).
        wavelength_data : numpy array
            wavelength data vector for each order. Shape (N_orders, N_pixels).
        fiber : int
            Fiber number.
        pm : PeakModeller
            Fit model for peaks.  If given, will be used while fitting lines, otherwise a Gaussian Model will be used.
        filename : str
            Filename of raw file.  Used for book-keeping.
        """
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
        """
        Returns the data.

        Returns:
            data (dataframe): Data.
        """
        return self._data

    @property
    def orders(self):
        """
        Returns the orders.

        Returns:
            orders (list): Orders.
        """
        return self._orders

    def normalize_orders(self, order, min_order, max_order):
        '''
        Converts the physical order to a normalized order.

        Parameters
        ----------
        order : int
            Physical order.
        min_order : int
            Minimum physical order.
        max_order : int
            Maximum physical order.

        Returns
        -------
        norm_order : float
            Normalized order.
        '''
        norm_order = (order - min_order)/(max_order - min_order) * 2. - 1.
        return norm_order

    def normalize_pixel(self, pixel):
        '''
        Converts the physical pixel to a normalized pixel in the range [-1, 1].
        The normalization is done based on the number of pixels in box_data.

        Parameters
        ----------
        pixel : int
            Physical pixel.

        Returns
        -------
        norm_pixel : float
            Normalized pixel.
        '''
        # norm_pixel = (pixel - min_pixel)/(max_pixel - min_pixel) * 2. - 1.
        norm_pixel = pixel / self.box_data.shape[1] * 2. - 1.
        return norm_pixel

    @deprecated(since="DRAGONS integration", message=_DEPRECATION_MSG)
    def data_flattened(self, box_data=False):
        '''
        Returns the data flattened.

        Parameters
        ----------
        box_data : bool
            If False, data is optimal extraction, if True, data is box extraction.
        Returns
        -------
        data : Tuple: (intensity, wavelength) as a single array
        '''
        data_selection = 'box_data' if box_data else 'opt_data'
        return np.hstack(self.data[data_selection].values), np.hstack(self.data['wavelength'].values)

    @deprecated(since="DRAGONS integration", message=_DEPRECATION_MSG)
    def min_wavelength(self):
        """
        Returns: Minimum wavelength in nm.
        """
        return np.min(self.data['wavelength'])

    @deprecated(since="DRAGONS integration", message=_DEPRECATION_MSG)
    def max_wavelength(self):
        """
        Returns: Maximum wavelength in nm.
        """
        return np.max(self.data['wavelength'])

    @deprecated(since="DRAGONS integration", message=_DEPRECATION_MSG)
    def find_orders_containing_wavelength(self, wavelength):
        """
        Finds the orders that contain the specified wavelength.
        TODO: Optimize this function.
        Args:
            wavelength (float): Wavelength to find orders for.

        Returns:
            found (list): List of orders that contain the specified wavelength.
        """
        found = []
        for i in self.data['wavelength']:
            if min(i) < wavelength < max(i):
                found.append(i)
        return found

    @deprecated(since="DRAGONS integration", message=_DEPRECATION_MSG)
    def blaze_correct(self, flat_spectrum, box_data = False):
        """
        Adds deblazed values in the Pandas dataframe for either box or optimal extraction.
        Args:
            flat_spectrum (MaroonXSpectrum object): Flat spectrum.
            box_data (bool): If false, data is optimal extraction, if true, data is box extraction

        Returns:
            None
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
        Applies the specified wavelength solution to the spectrum.

        Args:
            wavelength_solution (WavelengthSolution): Wavelength solution to apply.
        """
        log = self.log
        log.info('Applying wavelength solution')

        for o in self.orders:
            x = np.arange(len(self.data.loc[o]['box_data']))
            self._data.loc[o]["wavelength"] = wavelength_solution.get_wavelength(x, o)
