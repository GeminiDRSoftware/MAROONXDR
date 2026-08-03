"""
Spectrum class for fitting blaze functions.

``FlatSpectrum`` extends ``EchelleSpectrum`` with per-order spline
fitting of the blaze function of a box extracted master flat, which is
required for correcting echelle spectra.
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from time import time
from scipy.interpolate import LSQUnivariateSpline
from .echellespectrum import EchelleSpectrum


from astropy.utils.decorators import deprecated

_DEPRECATION_MSG = (
    "This method is inherited from the legacy MaroonX pipeline and "
    "is not used by the DRAGONS reduction; it may be removed in a "
    "future release."
)


class FlatSpectrum(EchelleSpectrum):
    r"""
    Extracted 1-D spectrum of a flat fiber.

    Extends ``EchelleSpectrum`` with per-order blaze function fitting.
    The blaze dictionaries, keyed by physical order, are empty until
    ``fit_blaze`` is called.

    Parameters
    ----------
    \*\*kwargs
        Keyword arguments passed on to the ``EchelleSpectrum``
        constructor.

    Attributes
    ----------
    blaze : dict
        Normalized (max = 1) blaze function per order, from a spline fit
        to the box extracted master flatfield.

    blaze_correction : dict
        Fitted spline object (``scipy.interpolate.LSQUnivariateSpline``)
        per order.

    blaze_norm : dict
        Normalization factor for the blaze per order (peak flux in each
        flatfield stripe).
    """
    def __init__ (self, **kwargs):
        super(FlatSpectrum, self).__init__(**kwargs)
        self.blaze = {}
        self.blaze_correction = {}
        self.blaze_norm = {}

    def fit_blaze(self, debug = 0, spline_kwargs = None, n_knots=50):
        """
        Fit the blaze function of each order and return it.

        Each order's box extracted data is fit with a least-squares
        spline, using zero weights for NaN and zero-valued pixels. After
        an outlier rejection pass (pixels deviating from the fit by more
        than 5% of the fit value are rejected) the spline is fit again
        and normalized to a maximum of 1. The results are stored in the
        ``blaze``, ``blaze_correction``, and ``blaze_norm`` dictionaries.
        Orders without valid data are skipped with a warning.

        Parameters
        ----------
        debug : int
            Debug level. Above 1, spline fit timings are logged; above 2,
            the fit of each order is plotted.

        spline_kwargs : dict
            Additional keyword arguments passed to
            ``scipy.interpolate.LSQUnivariateSpline``.

        n_knots : int
            Number of knots to use in the spline fit. Default is 50.

        Returns
        -------
        dict
            Normalized blaze function per order (the ``blaze`` attribute).
        """
        log = self.log
        if spline_kwargs is None:
            spline_kwargs = {}
        for o in self.orders:
            try:
                xx = np.arange(len(self.data.loc[o]['box_data']))
                t = np.rint(np.linspace(0, len(xx), n_knots))[1:-1]
                intensity = np.asarray(self.data.loc[o]['box_data'], dtype=float)
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
                plt.plot(self.data.loc[o]["box_data"], label='original data')
                intensity[intensity == 0] = np.nan
                plt.plot(intensity, label='data with outliers masked')
                plt.plot(self.blaze[o]*self.blaze_norm[o], label='Blaze fit')
                plt.show()
        return self.blaze

    @deprecated(since="DRAGONS-integration", message=_DEPRECATION_MSG)
    def save_blaze_function(self, filename):
        """
        Save the blaze function to a file as an ascii table.

        Parameters
        ----------
        filename : str
            Name of the file to save the blaze function to.
        """
        # Convert blaze function to table
        blaze = self.blaze
        blaze_table = []
        for o in blaze:
            blaze_table.append([o, blaze[o]])

        # Save table as astropy table
        blaze_table = ap.table.Table(blaze_table, names=['order', 'blaze'])
        blaze_table.write(filename, format='ascii', overwrite=True)



