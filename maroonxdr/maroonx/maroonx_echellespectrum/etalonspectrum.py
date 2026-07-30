"""
Spectrum class for etalon calibration fibers.

``EtalonSpectrum`` extends ``EchelleSpectrum`` with the fitted etalon
peak data and a Fabry-Perot etalon model that converts order numbers
to wavelengths. It is used by the wavelength solution primitives to
assign wavelengths and order numbers to the etalon peaks.
"""
import numpy as np

from lmfit import parameter
from scipy.interpolate import UnivariateSpline, BSpline
from scipy.signal import medfilt

from .echellespectrum import EchelleSpectrum

from astropy.utils.decorators import deprecated

_DEPRECATION_MSG = (
    "This method is inherited from the legacy MaroonX pipeline and "
    "is not used by the DRAGONS reduction; it may be removed in a "
    "future release."
)

c = 3e8
class EtalonSpectrum(EchelleSpectrum):
    r"""
    Extracted 1-D spectrum of an etalon calibration fiber.

    Extends ``EchelleSpectrum`` with the fitted etalon peak data and the
    etalon model parameters used to convert interference order numbers
    to wavelengths.

    Parameters
    ----------
    peak_data : DataFrame
        Fitted etalon peak parameters of this fiber, indexed by
        (order, center) with the index values also kept as columns.
        The methods of this class add the columns
        ``WAVELENGTH_BY_THAR``, ``DISPERSION_MPS``, ``M``,
        ``M_FRACTION``, and ``WAVELENGTH``.

    poly_data : DataFrame or None
        Fitted etalon polynomial data. Currently never passed by
        ``MXSpectrum`` and not used.

    etalon_peaks_symmetric : bool
        Whether the peaks were fit with equal left and right sigmas.
        Stored but currently unused; the symmetric peak fits are
        selected at the ``MXSpectrum`` level instead. Default is False.

    \*\*kwargs
        Keyword arguments passed on to the ``EchelleSpectrum``
        constructor.

    Attributes
    ----------
    etalon_parameters : lmfit.parameter.Parameters
        Etalon model parameters, initialized with a cavity thickness of
        9.9985 mm. The wavelength solution primitives overwrite them
        with the fitted parameters from the reference wavelength file.
    """
    def __init__(self, peak_data, poly_data=None,  etalon_peaks_symmetric = False, **kwargs):
        super().__init__(**kwargs)
        self.peak_data = peak_data
        self.poly_data = poly_data
        self.etalon_peaks_symmetric = etalon_peaks_symmetric
        self.etalon_parameters = self.generate_etalon_parameters(l=9.9985)


    def generate_etalon_parameters(self, l=10.001, n=1, theta=0):
        """
        Generate the etalon model parameters.

        Returns an lmfit parameter set describing the Fabry-Perot etalon
        equation. Only ``l`` is allowed to vary (between 9.9 and 10.1 mm)
        when the parameters are fit.

        Parameters
        ----------
        l : float
            Etalon thickness in mm.

        n : float
            Refractive index of the (vacuum) etalon gap. Default is 1.

        theta : float
            Angle of incidence in radians. Default is 0.

        Returns
        -------
        lmfit.parameter.Parameters
            The etalon parameters ``l``, ``n``, and ``theta``.
        """
        p = parameter.Parameters()
        p.add("l", l, vary=True, min=9.9, max=10.1)
        p.add("n", n, vary=False)
        p.add('theta', theta, vary=False)
        return p

    def make_b_spline_from_pars(self, kind=5):
        """
        Create the etalon dispersion B-spline from the etalon parameters.

        Collects the ``knot_*`` (knots) and ``disp_*`` (coefficients)
        entries of ``etalon_parameters`` into a B-spline.

        Parameters
        ----------
        kind : int
            Degree of the spline. Default is 5.

        Returns
        -------
        scipy.interpolate.BSpline
            Dispersion correction spline. Does not extrapolate outside
            the knot range.
        """
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
        """
        Convert order numbers to wavelengths with dispersion correction.

        Parameters
        ----------
        mm : int or ndarray
            Etalon interference order number(s) of the peak(s).

        Returns
        -------
        float or ndarray
            Wavelength in nm.
        """
        spl = self.make_b_spline_from_pars()
        parameters = self.etalon_parameters
        return (2. * (parameters['l'] - spl(1 / mm)*parameters['l']) * np.cos(parameters['theta']) * parameters['n'] / mm) * 1e6

    def peak_to_wavelength(self, mm):
        """
        Convert etalon peak order numbers to wavelengths.

        Uses the dispersion-corrected etalon equation when the
        dispersion spline parameters are present in
        ``etalon_parameters``, and the plain etalon equation otherwise.

        Parameters
        ----------
        mm : int or ndarray
            Etalon interference order number(s) of the peak(s).

        Returns
        -------
        float or ndarray
            Wavelength in nm.
        """
        parameters = self.etalon_parameters
        if 'knot_0' in parameters:
            return self.peak_to_wavelength_spline(mm)
        else:
            return (2.0 * (parameters["l"]) * np.cos(parameters["theta"]) * parameters["n"] / mm) * 1e6

    def peak_to_wavelength_nodispersion(self, mm):
        """
        Convert order numbers to wavelengths without dispersion correction.

        Parameters
        ----------
        mm : int or ndarray
            Etalon interference order number(s) of the peak(s).

        Returns
        -------
        float or ndarray
            Wavelength in nm.
        """
        parameters = self.etalon_parameters
        return (2.0 * (parameters["l"]) * np.cos(parameters["theta"]) * parameters["n"] / mm) * 1e6

    def guess_m(self, wl):
        """
        Guess the interference order numbers for the given wavelengths.

        Inverts the etalon equation and rounds to the nearest integer
        order number. When the dispersion spline parameters are present
        in ``etalon_parameters``, a first pass without dispersion
        correction provides the order numbers at which the dispersion
        correction is then evaluated.

        Parameters
        ----------
        wl : float or ndarray
            Wavelength(s) in nm.

        Returns
        -------
        tuple of ndarray
            Rounded interference order number(s) and the fractional
            remainder (exact minus rounded value).
        """
        parameters = self.etalon_parameters
        if 'knot_0' in parameters:
            spl = self.make_b_spline_from_pars()
            m0 = np.array(np.rint(2.0 * parameters['l'] * np.cos(parameters['theta']) *  parameters['n'] / (wl / 1e6)), dtype=int)
            m_float  = np.array(2.0 * (parameters['l']- spl(1/m0)*parameters['l']) * \
                                np.cos(parameters['theta']) *  parameters['n'] / (wl / 1e6))
            m_int = np.rint(m_float).astype(int)
            return m_int, m_float-m_int
        else:
            m_float = np.array(2.0 * parameters['l'] * np.cos(parameters['theta']) * parameters['n'] / (wl / 1e6))
            m_int = np.rint(m_float).astype(int)
            return m_int, m_float-m_int

    @deprecated(since="DRAGONS-integration", message=_DEPRECATION_MSG)
    def get_peak_data(self, order, data="all"):
        """
        Get the peak data for the specified order.

        Parameters
        ----------
        order : int
            Order to get the peak data for.

        data : str
            Column of ``peak_data`` to get, or "all" for the whole row
            selection.

        Returns
        -------
        DataFrame or ndarray
            Peak data of the specified order: a DataFrame for "all",
            otherwise the values of the selected column.
        """
        if self.peak_data is not None:
            if data == "all":
                return self.peak_data.loc[order, :]
            else:
                return self.peak_data.loc[order, data].values

    def apply_wavelength_solution(self, wavelength_solution):
        """
        Apply the specified wavelength solution to the spectrum and peaks.

        Delegates to ``EchelleSpectrum.apply_wavelength_solution`` for
        the per-order wavelength column, then additionally evaluates the
        solution at the etalon peak centers and stores the result in the
        ``WAVELENGTH_BY_THAR`` column of ``peak_data``. If the solution
        provides per-order mean residuals (``order_means``, in m/s),
        they are applied as a correction factor. Orders without peak data
        are skipped with a warning.

        Parameters
        ----------
        wavelength_solution : WavelengthSolution
            Wavelength solution to apply.
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
                peak_centers = self.peak_data.loc[order, "CENTER"].values
                # Get the wavelength for the peak centers
                wls = wavelength_solution.get_wavelength(peak_centers, order)
                self.peak_data.loc[order, "WAVELENGTH_BY_THAR"] = wls * (1.0 + correction/3e8)
            except KeyError:
                log.warning("No data for order {}".format(order))

    def apply_wavelength_vector(self, debug=0):
        """
        Assign wavelengths and dispersions to the etalon peaks.

        For each order, a cubic interpolating spline of the wavelength
        vector in ``data`` is evaluated at the peak centers and stored
        in the ``WAVELENGTH_BY_THAR`` column of ``peak_data``. The
        spline derivative at the centers gives the local dispersion,
        stored in m/s per pixel in the ``DISPERSION_MPS`` column, with
        the two edge values replaced by their nearest neighbour. Orders
        without peak data are skipped with a warning.

        Parameters
        ----------
        debug : int
            Debug level. Currently unused.
        """
        # Calculate the wavelength for peaks:
        for i, order in enumerate(self.orders):
            try:
                # y = self.wavelength_data[i]
                # keeps original (potentially unsorted) order
                y = np.asarray(self.data.loc[order, 'wavelength'])

                x = np.arange(len(y))
                if y[0] > y[-1]:
                    y = y[::-1]

                # Create a cubic spline for the wavelength data
                spl = UnivariateSpline(x, y, k=3, s=0)
                peak_centers = self.peak_data.loc[order, "CENTER"].values
                wls_vector = spl(peak_centers)

                # Calculate the dispersion in m/s
                splder = spl.derivative()
                dispersion = splder(peak_centers)
                dispersion[0] = dispersion[1]
                dispersion[-1] = dispersion[-2]
                dispersion_mps = dispersion * 3e8 / wls_vector

                # Assign wavelengths and dispersion to the peak data
                self.peak_data.loc[order, "WAVELENGTH_BY_THAR"] = wls_vector
                self.peak_data.loc[order, 'DISPERSION_MPS'] = dispersion_mps

            except KeyError:
                self.log.warning("No data for order {}".format(order))

    def guess_peak_numbers(self, debug=0, plot_title="", drop_outliers=True):
        """
        Guess the interference order number of each etalon peak.

        For each order, ``guess_m`` is applied to the
        ``WAVELENGTH_BY_THAR`` column of ``peak_data`` to fill the ``M``
        and ``M_FRACTION`` columns. Peaks whose residual between the
        ThAr-based and etalon-equation wavelengths deviates by more than
        5 sigma from a median-filtered baseline are logged and, if
        ``drop_outliers`` is True, removed. Jumps of up to four
        interference orders between neighbouring echelle orders are
        detected from the median residuals and corrected. Finally the
        ``WAVELENGTH`` column is filled with the etalon-equation
        wavelengths of the corrected order numbers.

        Parameters
        ----------
        debug : int
            Debug level. Currently unused.

        plot_title : str
            Title for the plot. Currently unused; plotting is done by
            ``plot_etalon_dispersion``.

        drop_outliers : bool
            If True (default), remove outlier peaks from ``peak_data``.
            Pass False when peak centers have been drift-corrected so
            that the MultiIndex still holds pre-correction center
            values: in that case the drop key (column ``CENTER``)
            differs from the index label (original ``CENTER``) and the
            drop would silently remove the wrong row. Legacy pandas
            1.0.1 raised KeyError in this situation but gets silently
            skipped via ``except: pass``, keeping all peaks.

        Returns
        -------
        DataFrame
            The updated ``peak_data``.
        """
        log = self.log
        # guess peak numbers based on wavelength and etalon model:
        for order in self.orders:

            # order_mask = self.peak_data["ORDER"] == order
            # if np.count_nonzero(order_mask) == 0:
            #     log.warning(f"No data for order {order}")
            #     continue

            self.peak_data.loc[order, 'M'],self.peak_data.loc[order, 'M_FRACTION'] = \
                self.guess_m(self.peak_data.loc[order, 'WAVELENGTH_BY_THAR'])

            # correct for inter-order jumps
            self.peak_data.sort_values("WAVELENGTH_BY_THAR", inplace=True, ascending=False)
            wl_peaks_by_thar = self.peak_data.loc[order, 'WAVELENGTH_BY_THAR'].values
            wl_by_etaloneq = self.peak_to_wavelength(self.peak_data.loc[order, 'M'].values)
            y = (wl_peaks_by_thar - wl_by_etaloneq) / wl_peaks_by_thar * c
            # m = self.peak_data.loc[order, 'M'].values
            residuals_flattened = y - medfilt(y, 11)
            bad = np.abs(residuals_flattened) > 5.0 * np.nanstd(residuals_flattened)

            if np.count_nonzero(bad) > 0:
                if np.count_nonzero(bad) < 5:
                    log.fullinfo(f'{np.count_nonzero(bad)} bad lines removed in order {order}')
                else:
                    log.warning(f'{np.count_nonzero(bad)} bad lines removed in order {order}')

                if drop_outliers:
                    dropindex = self.peak_data.loc[order].iloc[np.where(bad)].index
                    log.warning(f"Dropping {list(dropindex)} indices in order {order}")
                    for idx in dropindex:
                        self.peak_data.drop(index=(order, idx), inplace=True)

        # correct jumps between orders
        for order in (self.orders[1:])[::-1]:
            # order_mask = self.peak_data["ORDER"] == order
            # order_m1_mask = self.peak_data["ORDER"] == order - 1

            wl_peaks_by_thar = self.peak_data.loc[order-1, 'WAVELENGTH_BY_THAR'].values
            wl_by_etaloneq = self.peak_to_wavelength(self.peak_data.loc[order-1, 'M'].values)

            y2 = np.median((wl_peaks_by_thar - wl_by_etaloneq) / wl_peaks_by_thar * c)

            wl_peaks_by_thar = self.peak_data.loc[order, 'WAVELENGTH_BY_THAR'].values
            wl_by_etaloneq = self.peak_to_wavelength(self.peak_data.loc[order, 'M'].values)
            y1 = np.median((wl_peaks_by_thar - wl_by_etaloneq) / wl_peaks_by_thar * c)

            oldvalues = self.peak_data.loc[order-1, 'M'].values

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
                self.peak_data.loc[order-1, 'M'] = oldvalues + jump

            for order in self.orders:
                self.peak_data.loc[order, 'WAVELENGTH'] = self.peak_to_wavelength(self.peak_data.loc[order, 'M'].values)

        return self.peak_data

    def plot_etalon_dispersion(self, plot_title = "", plot_mfraction = True):
        """
        Plot the etalon dispersion and residuals.

        Parameters
        ----------
        plot_title : str
            Title for the plot.

        plot_mfraction : bool
            If True, the third panel shows the deviation from integer
            peak numbers; otherwise the residuals against normalized x.

        Returns
        -------
        matplotlib.figure.Figure
            The created figure.
        """
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(8, 8))
        fig.subplots_adjust(bottom=.07, left=0.14, right=0.96, top=0.95, hspace=0.40)
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312, sharex=ax1)
        ax3 = fig.add_subplot(313)

        ax2.set_title('Etalon residuals after dispersion correction ' + plot_title)
        ax2.set_xlabel('Wavelength (nm)')
        ax2.set_ylabel('Residuals (m/s)')
        residuals = (self.peak_to_wavelength \
                    (self.peak_data.loc[:, 'M'].values) - \
                    self.peak_data.loc[:, 'WAVELENGTH_BY_THAR'])
        residuals = residuals * c / self.peak_data.loc[:, 'WAVELENGTH_BY_THAR'] # Convert to m/s and normalize
        ax2.scatter(self.peak_data.loc[:, 'WAVELENGTH_BY_THAR'], residuals,
                    c=self.peak_data.loc[:, 'ORDER'], cmap='nipy_spectral',
                    rasterized=True, marker='.', s=2)
        ax2.set_ylim(np.nanmean(residuals) - 30,np.nanmean(residuals) + 30)
        ax2.text(0.6, 0.05,
                 f'std: {np.std(residuals):.1f} m/s, mean: {np.mean(residuals):.2f} m/s',
                 transform=ax2.transAxes)
        if 'knot_0' in self.etalon_parameters:
            dispersion = (self.peak_to_wavelength_nodispersion(self.peak_data.loc[:, 'M'].values) -
                          self.peak_to_wavelength(self.peak_data.loc[:, 'M'].values))
            dispersion = dispersion / self.peak_data.loc[:, 'WAVELENGTH_BY_THAR'] * c
            ax1.set_title('Etalon residuals ' + plot_title)
            ax1.set_xlabel('Wavelength (nm)')
            ax1.set_ylabel('Residuals (m/s)')
            uresiduals = (self.peak_to_wavelength_nodispersion(self.peak_data.loc[:, 'M'].values) -
                         self.peak_data.loc[:, 'WAVELENGTH_BY_THAR'])
            uresiduals = uresiduals * 3e8 / self.peak_data.loc[:,'WAVELENGTH_BY_THAR']
            ax1.scatter(self.peak_data.loc[:, 'WAVELENGTH_BY_THAR'], uresiduals,
                        c=self.peak_data.loc[:, 'ORDER'], cmap='nipy_spectral', rasterized=True, marker='.', s=2)
            ax1.plot(self.peak_data.loc[:, 'WAVELENGTH_BY_THAR'], dispersion, 'r-', rasterized=True)
            ax1.text(0.6, 0.05,
                     f'std: {np.std(uresiduals):.1f} m/s, mean: {np.mean(uresiduals):.2f} m/s',
                     transform=ax1.transAxes)

            if plot_mfraction:

                ax3.set_title('Deviation from integer peak # ' + plot_title)
                ax3.set_xlabel('Wavelength (nm)')
                ax3.set_ylabel('Fractional order #')
                ax3.set_xlim(ax2.get_xlim()[0], ax2.get_xlim()[1])
                ax3.set_ylim(ax2.get_ylim()[0] / 8000, ax2.get_ylim()[1] / 8000)
                ax3.scatter(self.peak_data.loc[:, 'WAVELENGTH_BY_THAR'], self.peak_data.loc[:, 'M_FRACTION'],
                            c=self.peak_data.loc[:, 'ORDER'], cmap='nipy_spectral', rasterized=True, marker='.', s=2)

            else:
                ax3.set_title('Etalon residuals after dispersion correction ' + plot_title)
                ax3.set_xlabel('X (normalized)')
                ax3.set_ylabel('Residuals (m/s)')
                ax3.set_xlim(-1,1)
                ax3.set_ylim(ax2.get_ylim()[0], ax2.get_ylim()[1])
                x = self.normalize_pixel(self.peak_data.loc[:, 'CENTER'])
                ax3.scatter(x, residuals,
                            c=self.peak_data.loc[:, 'ORDER'], cmap='nipy_spectral', rasterized=True, marker='.', s=2)

        return fig