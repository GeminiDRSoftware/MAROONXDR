
import numpy as np
import astropy.units as u
import math
from astroquery.nist import Nist
from numba import jit
from numpy.lib.stride_tricks import as_strided
import time as time
from gempy.utils import logutils

class Spectrum_Utils:
    def __init__(self):
        self.log = logutils.get_logger(__name__)

    def cross_correlation(self, x, y, maxlag):
        '''
        Cross correlation with a maximum number of lags.
        'x' and 'y' must be one-dimensional numpy arrays with the same length.

        This computes the same result as numpy.correlate(x, y, mode='full')[len(a)-maxlag-1:len(a)+maxlag].

        The return value has length 2*maxlag+1.

        Parameters
        ----------
        x : numpy array
            First array to correlate.
        y : numpy array
            Second array to correlate.
        maxlag : int
            Maximum number of lags.

        Returns
        -------
        xcorr : numpy array
            Cross correlation.
        '''
        log = self.log
        if len(x) != len(y):
            log.error('x and y must have the same length.')
            return
        pixel_y = np.pad(y.conj(), 2*maxlag, mode='constant')
        T = as_strided(pixel_y[maxlag:], shape=(2*maxlag+1, len(y)+2*maxlag), strides=(-pixel_y.strides[0], pixel_y.strides[0]))
        pixel_x = np.pad(x, maxlag, mode='constant')
        return T.dot(pixel_x)

    @jit(nopython=True, cache=True)
    def remove_jumps(self, interference_number, residuals):
        '''
        Removes jumps in the residuals.

        Parameters
        ----------
        interference_number : numpy array
            Interference number.
        residuals : numpy array
            Residuals.
        '''
        log = self.log
        i_num_range = np.arange(len(interference_number))
        # Number of jumps
        num_jumps = 0
        for i in range(len(residuals) - 1):
            if residuals[i+1] - residuals[i] > 40000:
                interference_number[i_num_range > i] -= 4
                num_jumps += 1
            elif residuals[i + 1] - residuals[i] > 25000:
                interference_number[i_num_range > i] -= 3
                num_jumps += 1
            elif residuals[i + 1] - residuals[i] > 15000:
                interference_number[i_num_range > i] -= 2
                num_jumps += 1
            elif residuals[i + 1] - residuals[i] > 3000:
                interference_number[i_num_range > i] -= 1
                num_jumps += 1

            if residuals[i + 1] - residuals[i] < -40000:
                interference_number[i_num_range > i] += 4
                num_jumps += 1
            elif residuals[i + 1] - residuals[i] < -25000:
                interference_number[i_num_range > i] += 3
                num_jumps += 1
            elif residuals[i + 1] - residuals[i] < -15000:
                interference_number[i_num_range > i] += 2
                num_jumps += 1
            elif residuals[i + 1] - residuals[i] < -3000:
                interference_number[i_num_range > i] += 1
                num_jumps += 1

            log.debug('Number of jumps: {}'.format(num_jumps))
            return interference_number, num_jumps

    def pull_catalogue_lines(min_wl, max_wl, catalogue = 'Th', wavelength_type = 'vacuum'):
        """
        Pulls lines from the NIST catalogue in the specified wavelength range.

        Parameters
        ----------
        min_wl : float
            Minimum wavelength.
        max_wl : float
            Maximum wavelength.
        catalogue : str
            Catalogue to pull lines from.  Can be 'Th' or 'I2'.
        wavelength_type : str
            Wavelength type.  Can be 'vacuum' or 'air'.

        Returns
        -------
        lines : numpy array
            Wavelengths of the lines.
        """
        table_lines = Nist.query(min_wl * u.nm, max_wl * u.nm, linename = catalogue, output_order = 'wavelength',
                             wavelength_type = wavelength_type)

        intensity_lines = [table_lines('Rel.')]
        error_count = 0
        for pos, val in enumerate(intensity_lines[0]):
            try:
                new_val = float(val)
            except ValueError:
                new_val = 0.0
                error_count += 1
            if not (math.isnan(new_val) or math.isinf(new_val)):
                intensity_lines[0][pos] = new_val

        intensity_lines = np.array(intensity_lines[0], dtype=float)
        th_index = intensity_lines > 0.0
        th = np.array(table_lines('Ritz'))
        return th[th_index], intensity_lines[th_index]