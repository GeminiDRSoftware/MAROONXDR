import numpy as np
import matplotlib.pyplot as plt
from .echellespectrum import EchelleSpectrum
from lmfit import parameter
from scipy.interpolate import UnivariateSpline, BSpline
from scipy.signal import medfilt

c = 3e8
class EtalonSpectrum(EchelleSpectrum):
    """
    Describes an etalon spectrum.  Inherits from EchelleSpectrum.
    """
<<<<<<< HEAD
    def __init__(self, peak_data, poly_data,  etalon_peaks_symmetric = False, **kwargs):
=======
    def __init__(self, peak_data, **kwargs):
>>>>>>> 9289dd70091d33872455b5dab3b5248b5a783cbe
        """
        Initializes the EtalonSpectrum object.
        """
        super().__init__(**kwargs)
        self.peak_data = peak_data
<<<<<<< HEAD
        self.poly_data = poly_data
        self.etalon_peaks_symmetric = etalon_peaks_symmetric
=======
>>>>>>> 9289dd70091d33872455b5dab3b5248b5a783cbe
        self.etalon_parameters = self.generate_etalon_parameters(l=9.9985)


    def generate_etalon_parameters(self, l=10.001, n=1, theta=0):
        """
        Generates the etalon parameters.
        Args:
            l (float): Etalon thickness in mm.
            n (int): Etalon order.
            theta (float): Etalon angle in degrees.
        """
        p = parameter.Parameters()
        p.add("l", l, vary=True, min=9.9, max=10.1)
        p.add("n", n, vary=False)
        p.add('theta', theta, vary=False)
        return p

    def make_b_spline_from_pars(self, kind=5):
        '''
        Creates a B-spline from the parameters.

        Args:
            kind (int): Order of the spline.

        Returns:
            spline (scipy.interpolate.BSpline): B-spline.
        '''
        disp_params = []
        for par in self.etalon_parameters:
            if "disp_" in par:
                disp_params.append(self.etalon_parameters[par].value)
        disp_params = np.array(disp_params, dtype=float)
        t = []
        for par in self.etalon_parameters:
            if "knot_" in par:
                t.append(self.etalon_parameters[par].value)
        knots = np.array(t, dtype=float)

        return BSpline(knots, disp_params, kind, extrapolate=False)

    def peak_to_wavelength_spline(self, mm):
        '''
        Creates a spline to convert etalon peaks to wavelengths

        Args:
            mm (float): Etalon peak in mm.

        Returns:
            wavelength (float): Wavelength in nm.
        '''
        spl = self.make_b_spline_from_pars()
        parameters = self.etalon_parameters
        return (2. * (parameters['l'] - spl(1 / mm)*parameters['l']) * np.cos(parameters['theta']) * parameters['n'] / mm) * 1e6

    def peak_to_wavelength(self, mm):
        '''
        Applies the spline to convert etalon peaks to wavelengths

        Args:
            mm (float): Etalon peak in mm.
        Returns:
            wavelength (float): Wavelength in nm.
        '''
        parameters = self.etalon_parameters
        if 'knot_0' in parameters:
            return self.peak_to_wavelength_spline(mm)
        else:
            return (2.0 * (parameters["l"]) * np.cos(["theta"]) * parameters["n"] / mm) * 1e6

    def peak_to_wavelength_nodispersion(self, mm):
        '''
        Applies the spline to convert etalon peaks to wavelengths, but without dispersion correction

        Args:
            mm (float): Etalon peak in mm.
        Returns:
            wavelength (float): Wavelength in nm.
        '''
        parameters = self.etalon_parameters
        return (2.0 * (parameters["l"]) * np.cos(parameters["theta"]) * parameters["n"] / mm) * 1e6

    def guess_m(self, wl):
        '''
        Use the etalon parameters to guess the order number for a given wavelength.
        '''
        parameters = self.etalon_parameters
        if 'knot_0' in parameters:
            spl = self.make_b_spline_from_pars()
            m0 = np.array(np.rint(2.0 * parameters['l'] * np.cos(parameters['theta']) *  parameters['n'] / (wl / 1e6)), dtype=int)
            m_float  = np.array(2.0 * (parameters['l']- spl(1/m0)*parameters['l']) * \
                                np.cos(parameters['theta']) *  parameters['n'] / (wl / 1e6))
            m_int = np.rint(m_float).astype(np.int)
            return m_int, m_float-m_int
        else:
            m_float = np.array(2.0 * parameters['l'] * np.cos(parameters['theta']) * parameters['n'] / (wl / 1e6))
            m_int = np.rint(m_float).astype(np.int)
            return m_int, m_float-m_int

    def get_peak_data(self, order, data="all"):
        """
        Gets the peak data for the specified order.

        Args:
            order (int): Order to get the peak data for.
            data (str): Type of data to get.  Can be "all", "box", or "optimal".

        Returns:
            peak_data (1D array) : Peak data for the specified order.
        """
        if self.peak_data is not None:
            if data == "all":
                return self.peak_data.loc[order, :]
            else:
                return self.peak_data.loc[order, data].values

    def apply_wavelength_solution(self, wavelength_solution):
        """
        Applies the specified wavelength solution to the spectrum.

        Args:
            wavelength_solution (WavelengthSolution): Wavelength solution to apply.
        """
        super().apply_wavelength_solution(wavelength_solution)
        log = self.log
        # Calculate the wavelength for peaks:
        for order in self.orders:
            try:
                if wavelength_solution.order_means is not None:
                    correction = wavelength_solution.order_means[order]
                else:
                    correction = 0
                self.peak_data.loc[order, "wavelength_by_thar"] = wavelength_solution.get_wavelength(
                    self.peak_data.loc[order, "center"], order) * (1.0 + correction/3e8)
            except KeyError:
                log.warning("No data for order {}".format(order))

    def apply_wavelength_vector(self, debug = 0):
        """
        Assign wavelengths to each Etalon peak based on a cubic spline interpolation of the wavelength vector
        in the data.

        Args:
            debug (int): Debug level.

        Returns:
            None
        """
        # Calculate the wavelength for peaks:
        for order in self.orders:
            try:
                x = np.arange(len(self.data.loc[order]['wavelength']))
                y = self.data.loc[order]['wavelength']
                if y[0] > y[-1]:
                    y = y[::-1]
                spl = UnivariateSpline(x, y, k=3, s=0)
                self.peak_data.loc[order, "wavelength_by_spline"] = spl(self.peak_data.loc[order, "center"])
            except KeyError:
                self.log.warning("No data for order {}".format(order))

    def guess_peak_numbers(self, debug=0, plot_title = ""):
        """
        Guess peak numbers based on the wavelength and etalon model.

        Args:
            debug (int): Debug level.
            plot_title (str): Title for the plot.

        Returns:
            peak_numbers (1D array): Peak numbers.
        """
        log = self.log
        # guess peak numbers based on wavelength and etalon model:
        for order in self.orders:
            try:
                self.peak_data.loc[order, 'm'],self.peak_data.loc[order, 'm_fraction'] = \
                    self.guess_m(self.peak_data.loc[order, 'wavelength_by_thar'])

                # correct for inter-order jumps
                self.peak_data.sort_values("wavelength_by_thar", inplace=True, ascending=False)
                wl_peaks_by_thar = self.peak_data.loc[order, 'wavelength_by_thar'].values
                wl_by_etaloneq = self.peak_to_wavelength(self.peak_data.loc[order, 'm'].values)
                y = (wl_peaks_by_thar - wl_by_etaloneq) / wl_peaks_by_thar * c
                m = self.peak_data.loc[order, 'm'].values
                residuals_flattened = y - medfilt(y, 11)
                bad = np.abs(residuals_flattened) > 5.0 * np.nanstd(residuals_flattened)
                if np.count_nonzero(bad) > 0:
                    if np.count_nonzero(bad) < 5:
                        log.info(f'{np.count_nonzero(bad)} bad lines removed in order {order}')
                    else:
                        log.warning(f'{np.count_nonzero(bad)} bad lines removed in order {order}')
                    dropindex = list((self.peak_data.loc[order].iloc[np.where(bad)].center.values))
                    for do in dropindex:
                        self.peak_data.drop(index=(order,do),inplace=True)
            except:
                log.warning(f"No data for order {order}")

        # correct jumps between orders
        for order in (self.orders[1:])[::-1]:
            wl_peaks_by_thar = self.peak_data.loc[order-1, 'wavelength_by_thar'].values
            wl_by_etaloneq = self.peak_to_wavelength(self.peak_data.loc[order-1, 'm'].values)

            y2 = np.median((wl_peaks_by_thar - wl_by_etaloneq) / wl_peaks_by_thar * c)

            wl_peaks_by_thar = self.peak_data.loc[order, 'wavelength_by_thar'].values
            wl_by_etaloneq = self.peak_to_wavelength(self.peak_data.loc[order, 'm'].values)
            y1 = np.median((wl_peaks_by_thar - wl_by_etaloneq) / wl_peaks_by_thar * c)

            oldvalues = self.peak_data.loc[order-1, 'm'].values

            jump = 0

            if y2-y1 > 40000:
                jump = -4
            elif y2-y1 > 25000:
                jump = -3
            elif y2-y1 > 15000:
                jump = -2
            elif y2-y1 > 3000:
                jump = -1
            if y2 - y1 < -40000:
                jump = +4
            elif y2 - y1 < -25000:
                jump = +3
            elif y2-y1 < -15000:
                jump = +2
            elif y2-y1 < -3000:
                jump = +1

            if jump !=0:
                log.warning(f'Jump by {jump} IOs found and corrected between order {order-1} and order {order}')
                self.peak_data.loc[order - 1, 'm'] = oldvalues + jump

            for order in self.orders:
                self.peak_data.loc[order, 'wavelength'] = self.peak_to_wavelength(self.peak_data.loc[order, 'm'].values)

        if debug > 0:
            fig = self.plot_etalon_dispersion(plot_title)
            plt.show() # Show the plot
        return self.peak_data

    def plot_etalon_dispersion(self, plot_title = "", plot_mfraction = True):
        """
        Plots the etalon dispersion.

        Args:
            plot_title (str): Title for the plot.

        Returns:
            None
        """
        fig = plt.figure(figsize=(8, 8))
        fig.subplots_adjust(bottom=.07, left=0.14, right=0.96, top=0.95, hspace=0.40)
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312, sharex=ax1)
        ax3 = fig.add_subplot(313)

        ax2.set_title('Etalon residuals after dispersion correction ' + plot_title)
        ax2.set_xlabel('Wavelength (nm)')
        ax2.set_ylabel('Residuals (m/s)')
        residuals = (self.peak_to_wavelength \
                    (self.peak_data.loc[:, 'm'].values) - \
                    self.peak_data.loc[:, 'wavelength_by_thar'])
        residuals = residuals * c / self.peak_data.loc[:, 'wavelength_by_thar'] # Convert to m/s and normalize
        ax2.scatter(self.peak_data.loc[:, 'wavelength_by_thar'], residuals,
                    c=self.peak_data.loc[:, 'order'], cmap='nipy_spectral',
                    rasterized=True, marker='.', s=2)
        ax2.text(0.6, 0.05,
                 f'std: {np.std(residuals):.1f} m/s, mean: {np.mean(residuals):.2f} m/s',
                 transform=ax2.transAxes)
        if 'knot_0' in self.etalon_parameters:
            dispersion = (self.peak_to_wavelength_nodispersion(self.peak_data.loc[:, 'm'].values) -
                          self.peak_to_wavelength(self.peak_data.loc[:, 'm'].values))
            dispersion = dispersion / self.peak_data.loc[:, 'wavelength_by_thar'] * c
            ax1.set_title('Etalon residuals ' + plot_title)
            ax1.set_xlabel('Wavelength (nm)')
            ax1.set_ylabel('Residuals (m/s)')
            uresiduals = (self.peak_to_wavelength_nodispersion(self.peak_data.loc[:, 'm'].values) -
                         self.peak_data.loc[:, 'wavelength_by_thar'])
            uresiduals = uresiduals * 3e8 / self.peak_data.loc[:,'wavelength_by_thar']
            ax1.scatter(self.peak_data.loc[:, 'wavelength_by_thar'], uresiduals,
                        c=self.peak_data.loc[:, 'order'], cmap='nipy_spectral', rasterized=True, marker='.', s=2)
            ax1.plot(self.peak_data.loc[:, 'wavelength_by_thar'], dispersion, 'r-', rasterized=True)
            ax1.text(0.6, 0.05,
                     f'std: {np.std(residuals):.1f} m/s, mean: {np.mean(residuals):.2f} m/s',
                     transform=ax1.transAxes)

            if plot_mfraction:

                ax3.set_title('Deviation from integer peak # ' + plot_title)
                ax3.set_xlabel('Wavelength (nm)')
                ax3.set_ylabel('Fractional order #')
                ax3.set_xlim(ax2.get_xlim()[0], ax2.get_xlim()[1])
                ax3.set_ylim(ax2.get_ylim()[0] / 8000, ax2.get_ylim()[1] / 8000)
                ax3.scatter(self.peak_data.loc[:, 'wavelength_by_thar'], self.peak_data.loc[:, 'm_fraction'],
                            c=self.peak_data.loc[:, 'order'], cmap='nipy_spectral', rasterized=True, marker='.', s=2)

            else:
                ax3.set_title('Etalon residuals after dispersion correction ' + plot_title)
                ax3.set_xlabel('X (normalized')
                ax3.set_ylabel('Residuals (m/s)')
                ax3.set_xlim(-1,1)
                ax3.set_ylim(ax2.get_ylim()[0], ax2.get_ylim()[1])
                x = self.normalize_x(self.peak_data.loc[:, 'center'])
                ax3.scatter(x, residuals,
                            c=self.peak_data.loc[:, 'order'], cmap='nipy_spectral', rasterized=True, marker='.', s=2)

        return fig