"""
FlatSpectrum class for handling fitting of blaze functions required for correcting echelle spectra.
Inherits from EchelleSpectrum.
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from time import time
from scipy.interpolate import LSQUnivariateSpline
from .echellespectrum import EchelleSpectrum


class FlatSpectrum(EchelleSpectrum):
    def __init__ (self, **kwargs):
        super(FlatSpectrum, self).__init__(**kwargs)
        self.blaze = {}
        self.blaze_correction = {}
        self.blaze_norm = {}

    def fit_blaze(self, debug = 0, spline_kwargs = None, n_knots=50):
        """
        Finds blaze function and returns it.

        Args:
            debug (int): Debug level.
            spline_kwargs (dict): Dictionary of kwargs to pass to the spline fit.
            n_knots (int): Number of knots to use in the spline fit.

        Returns:
            blaze (1D array): Blaze function.
        """
        log = self.log
        if n_knots is None:
            n_knots = 50
        if spline_kwargs is None:
            spline_kwargs = {}
        for o in self.orders:
            try:
                xx = np.arange(len(self.data.loc[o]['box_intensity']))
                t = np.rint(np.linspace(0, len(xx), n_knots))[1:-1]
                intensity = self.data.loc[o]['box_intensity'].copy()
                w = np.any((np.isnan(intensity), (intensity == 0)), axis=0)  # weights where the data are NAN or truly zero
                intensity[w] = 0
                t_start = time()
                self.blaze_correction[o] = LSQUnivariateSpline(xx, intensity, t, **spline_kwargs, check_finite=True, w=~w)
                t_end = time()
                if debug > 1:
                    log.debug(f"Blaze correction for order {o} took {t_end-t_start:.2f} seconds")
                self.blaze[o] = self.blaze_correction[o](xx)
                # simply outlier rejection (data-blazefit > 5% of blazefit)
                intensity[np.abs(intensity-self.blaze[o]) > (self.blaze[o]*0.05)] = 0
                w = (intensity == 0)
                t_start = time()
                self.blaze_correction[o] = LSQUnivariateSpline(xx, intensity, t, **spline_kwargs, check_finite=True, w=~w)
                t_end = time()
                if debug > 1:
                    log.debug(f"Blaze correction for order {o} after outlier rejection took {t_end-t_start:.2f} seconds")
                self.blaze[o] = self.blaze_correction[o](xx)
                self.blaze_norm[o] = np.max(self.blaze[o])
                self.blaze[o] = self.blaze[o]/self.blaze_norm[o]
            except TypeError:
                log.warning(f"Order {o} has no valid data. No blaze function is calculated")
            if debug > 2:
                plt.figure()
                plt.title(f'Order: {o}')
                plt.plot(self.data.loc[o]["box_intensity"], label='original data')
                intensity[intensity == 0] = np.nan
                plt.plot(intensity, label='data with outliers masked')
                plt.plot(self.blaze[o]*self.blaze_norm[o], label='Blaze fit')
                plt.show()
        return self.blaze

    def save_blaze_function(self, filename):
        """
        Saves the blaze function to a file.

        Args:
            filename (str): Name of the file to save the blaze function to.

        Returns:
            None. File is saved to disk.
        """
        # Convert blaze function to table
        blaze = self.blaze
        blaze_table = []
        for o in blaze:
            blaze_table.append([o, blaze[o]])

        # Save table as astropy table
        blaze_table = ap.table.Table(blaze_table, names=['order', 'blaze'])
        blaze_table.write(filename, format='ascii', overwrite=True)



