'''
This class contains the echelle spectrum, from which all other spectrum types inherit.
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .spectrum_utils import Spectrum_Utils
from scipy import interpolate
from lmfit.models import GaussianModel, LinearModel
from gempy.gemini import gemini_tools as gt


class EchelleSpectrum:
    '''
    The echelle spectrum class contains all information about extracted 1-D echelle spectra.
    Each object contains several orders, and corresponding wavelength and intensity data.
    '''
    def __init__(self, orders, box_data=None, box_error=None
                 , opt_data=None, opt_error=None, wavelength_data=None, fiber=1, pm=None,
                 filename = None):
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
<<<<<<< HEAD

        # physical orders - sorted like the index
        self.orders = orders
=======
        self.box_data = box_data
        self.box_error = box_error
        self.opt_data = opt_data
        self.opt_error = opt_error
        self.filename = filename

        #Construct a dataframe with the data
        self.data = pd.DataFrame.from_dict({'box_intensity': box_data, 'box_error': box_error,
                                            'opt_intensity': opt_data, 'opt_error': opt_error,
                                            'wavelength': wavelength_data})
        #'wavelength': wavelength_data
        self.data.set_index(orders, inplace=True)
        self.data.sort_index(inplace=True)

        # physical orders - sorted like the index
        self.orders = orders

>>>>>>> 9289dd70091d33872455b5dab3b5248b5a783cbe
        min_order = np.min(orders)
        max_order = np.max(orders)

        # normalized orders
        self.norm_orders = self.normalize_orders(orders, min_order, max_order)

        self.etalon_parameters = None
        self.model = None

        self.fiber = fiber
        self.pm = pm
        self.log = gt.logutils.get_logger(__name__)

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

    def normalize_pixel(self, pixel, min_pixel, max_pixel):
        '''
        Converts the physical pixel to a normalized pixel.

        Parameters
        ----------
        pixel : int
            Physical pixel.
        min_pixel : int
            Minimum physical pixel.
        max_pixel : int
            Maximum physical pixel.

        Returns
        -------
        norm_pixel : float
            Normalized pixel.
        '''

        norm_pixel = (pixel - min_pixel)/(max_pixel - min_pixel) * 2. - 1.
        return norm_pixel

<<<<<<< HEAD
    '''
=======
>>>>>>> 9289dd70091d33872455b5dab3b5248b5a783cbe
    @property
    def data(self):
        """
        Returns the data.

        Returns:
            data (dataframe): Data.
        """
        return self.data
<<<<<<< HEAD
    '''
=======

>>>>>>> 9289dd70091d33872455b5dab3b5248b5a783cbe
    @property
    def orders(self):
        """
        Returns the orders.

        Returns:
            orders (list): Orders.
        """
        return self.orders

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

    def min_wavelength(self):
        """
        Returns: Minimum wavelength in nm.
        """
        return np.min(self.data['wavelength'])

    def max_wavelength(self):
        """
        Returns: Maximum wavelength in nm.
        """
        return np.max(self.data['wavelength'])

    def plot_as_2d(self, offset_orders=True, box_data=False):
        """
        Plots all orders in one plot,, but offset, so they do not overplot each other.

        Args:
            offset_orders (bool): If true, plots will be offset so they don't overplot
            box_data (bool): If false, data is optimal extraction, if true, data is box extraction

        Returns:
            None
        """
        data_selection = 'box_data' if box_data else 'opt_data'
        color = list(plt.cm.tab20(np.linspace(0, 1, 20)))
        styles = ['-', '--', '-.', ':']
        plt.figure()
        max_int  = np.max(np.max(self.data.loc[self.orders[0]][data_selection]))

        for o in self.orders:
            if max_int < np.max(self.data.loc[o][data_selection]):
                max_int = np.max(self.data.loc[o][data_selection])

        for i, o in enumerate(self.orders):
            offset = 0
            if offset_orders:
                offset = o
            else:
                max_int = 1
            plt.plot(self.data.loc[o][data_selection] / max_int + offset, color=color[i % 20],
                    linestyle=styles[i % 4], label='Order {}'.format(o))

        plt.legend(bbox_to_anchor=(1.01, 0.5), loc="center left")

    def plot_orders(self, orders, box_data=False):
        """
        Plots the specified orders in one plot.

        Args:
            orders (list): List of orders to plot.
            box_data (bool): If false, data is optimal extraction, if true, data is box extraction

        Returns:
            None
        """
        data_selection = 'box_data' if box_data else 'opt_data'
        plt.figure()
        if isinstance(orders, int):
            orders = [orders]
        for o in orders:
            plt.plot(self.data.loc[o][data_selection], label='Order {}'.format(o))
        plt.show()

    def plot_overlap(self, orders, data_selection='opt_intensity'):
        """
        Plots the overlap between the specified orders.

        Args:
            orders (list): List of orders to plot.
            data_selection (str): Type of data to plot.  Can be 'opt_intensity' or 'box_intensity'.

        Returns:
            None
        """
        if orders is None:
            orders = self.orders

        orders = np.unique(np.array(orders))
        orders.sort()
        for i in range(len(orders)-1):
            plt.figure()
            wl1 = self.data.loc[orders[i]]['wavelength']
            wl2 = self.data.loc[orders[i+1]]['wavelength']
            idx1 = wl1 < np.max(wl2)
            idx2 = wl2 > np.min(wl1)
            data1 = self.data.loc[orders[i]][data_selection]
            data2 = self.data.loc[orders[i+1]][data_selection]

            plt.plot(wl1[idx1], data1[idx1], label=f"Order{orders[i]}")
            plt.plot(wl2[idx2], data2[idx2], label=f"Order{orders[i+1]}")
            plt.legend()
        plt.show()

    def plot_orders_vs_wavelength(self, orders=None, data_selection='box_intensity', plot_title=''):
        """
        Plots the specified orders vs wavelength.

        Args:
            orders (list): List of orders to plot.
            data_selection (str): Type of data to plot.  Can be 'opt_intensity' or 'box_extraction'.

        Returns:
            None
        """
        if orders is None:
            orders = self.orders
        if isinstance(orders, int):
            orders = [orders]

        fig = plt.figure(figsize=(8, 8))
        fig.subplots_adjust(bottom=.07, left=0.14, right=0.96, top=0.95, hspace=0.40)
        ax = fig.add_subplot(311)
        ax.set_title(f'ThAr Lines ({data_selection}) '+plot_title)
        ax.set_ylabel('Counts (DN)')
        ax.set_xlabel('Wavelength (nm)')
        for o in orders:
            plt.plot(self.data.loc[o]['wavelength'], self.data.loc[o][data_selection],rasterized=True)

        return fig

    def plot_xcorr_orderoverlap(self, orders=None, data_selection='box_extraction',
                                fig=None, plot_title=''):
        """
        Plots the cross-correlation between the specified orders.

        Args:
            orders (list): List of orders to plot.
            data_selection (str): Type of data to plot.  Can be 'opt_intensity' or 'box_extraction'.
            fig (Figure): Figure to plot on.
            plot_title (str): Title for the plot.

        Returns:
            None
        """
        offsets = []
        if orders is None:
            orders = self.orders

        saturation = 250000.0
        whs = 10

        if fig is None:
            fig = plt.figure()
            ax0 = plt.subplot(211)
            ax1 = plt.subplot(212)
        else:
            ax0 = plt.subplot(312)
            ax1 = plt.subplot(313)
        ax0.set_title('X-Correlation for order-overlap regions '+plot_title)
        ax0.set_xlabel('Order x (vs order x+1)')
        ax0.set_ylabel('FWHM (m/s)')
        ax1.set_title('X-corr for order-overlap regions ' + plot_title)
        ax1.set_xlabel('Order x (vs order x+1)')
        ax1.set_ylabel('Shift (m/s)')

        for o in orders[:-1]:
            data1 = self.data.loc[o][data_selection]
            data2 = self.data.loc[o+1][data_selection]
            wave1 = self.data.loc[o]['wavelength']
            wave2 = self.data.loc[o+1]['wavelength']
            data1[data1 == 0] = np.nan
            data1[data2 == 0] = np.nan

            orderlength = len(data1)

            saturated1 = np.asarray(np.where(data1 > saturation)[0])
            saturated2 = np.asarray(np.where(data2 > saturation)[0])

            if saturated1.size == 0:
                if saturated2.size != 0:
                    saturated1 = saturated2
            if saturated2.size == 0:
                if saturated1.size != 0:
                    saturated2 = saturated1

            saturated = np.unique(np.concatenate([saturated1, saturated2]))

            if saturated.size > 0:
                for x in saturated:
                    left = x - whs if (x - whs) >= 0 else 0
                    right = x + whs if (x + whs) <= orderlength - 1 else orderlength - 1
                    data1[left:right] = np.nan
                    data2[left:right] = np.nan

                wave1_trimmed = wave1[:np.argmax(wave1>np.max(wave2))]
                wave2_trimmed = wave2[np.argmax(wave2 > np.min(wave1)):]

                max_overlap_length = max(len(wave1_trimmed),len(wave2_trimmed))

                xlog = np.logspace(np.log10(np.max([wave1[0],wave2[0]])), np.log10(np.min([wave1[-1],wave2[-1]])), 5 * max_overlap_length)
                mpspp = (xlog[1] - xlog[0]) / xlog[0] * 3e8

                mask = np.isnan(data1)
                f = interpolate.interp1d(wave1[~mask], data1[~mask], kind='slinear', fill_value='extrapolate')
                intensity1_log = f(xlog)

                mask = np.isnan(data2)
                f = interpolate.interp1d(wave2[~mask], data2[~mask], kind='slinear', fill_value='extrapolate')
                intensity2_log = f(xlog)

                if o > 173 :
                    fig = plt.figure()
                    plt.plot(xlog, intensity1_log)
                    plt.plot(xlog, intensity2_log)
                    plt.title(o)

                lag = 7000 // mpspp
                xcorr = Spectrum_Utils.cross_correlation(intensity1_log, intensity2_log, lag)

                corrx = np.arange(2 * lag + 1) - lag

                pars = LinearModel().make_params(intercept=xcorr.min(), slope=1e5)
                pars += GaussianModel().guess(xcorr - np.min(xcorr), x=corrx)
                mod = LinearModel() + GaussianModel()
                out = mod.fit(xcorr, pars, x=corrx)
                if o > 173:
                    fig = plt.figure()
                    print(out.fit_report())
                    out.plot_fit()
                    plt.show(block=False)

                # Add the dispersion in m/s/per-pixel for the log-spaced oversampled wavelength vector
                # This will be used later to convert the lag into an RV shift in m/s
                out.params.add('MPSPP', value=mpspp, vary=False)

                offsets.append(out.params['center'] * mpspp)
                if out.params['fwhm'].stderr is not None:
                    ax0.errorbar(o, out.params['fwhm'].value * mpspp, yerr= out.params['fwhm'].stderr * mpspp)
                else:
                    ax0.errorbar(o, out.params['fwhm'].value * mpspp, yerr=0)

                ax0.scatter(o, out.params['fwhm'].value * mpspp, marker='+')

                if out.params['center'].stderr is not None:
                    ax1.errorbar(o, out.params['center'].value * mpspp, yerr= out.params['center'].stderr * mpspp)
                else:
                    ax1.errorbar(o, out.params['center'].value * mpspp, yerr=0)

                ax1.scatter(o, out.params['center'].value * mpspp, marker='+')

        ax1.text(0.05,0.85,f'Median: {np.median(offsets):.1f} m/s', transform=ax1.transAxes)

        return fig

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
            self.data["blaze_corrected_"+data_selection] = deblazed_int
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
            x = np.arange(len(self.data.loc[o]['box_intesity']))
            self.data.loc[o]["wavelength"] = wavelength_solution.get_wavelength(x, o)
