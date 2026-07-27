"""
Class that describes a wavelength solution for a given spectrum.
"""
import numpy as np
import time
import warnings

from astropy.modeling import models, fitting

from gempy.utils import logutils


from astropy.utils.decorators import deprecated

_DEPRECATION_MSG = (
    "This method is inherited from the legacy MaroonX pipeline and "
    "is not used by the DRAGONS reduction; it may be removed in a "
    "future release."
)


class WavelengthSolution:
    def __init__(self, x_norm, orders, weights, wavelengths, poly_deg_x=5, poly_deg_y=5, max_x=3984):
        '''
        Initializes the WavelengthSolution object.

        Parameters
        ----------
        x_norm : 1D array
            Normalized x values.
        orders : 1D array
            Vector of physical echelle orders. Same length as x and wavelengths.
        weights : 1D array
            Weights for each pixel.
        wavelengths : 1D array
            Wavelengths of lines
        poly_deg_x : int
            Degree of Legendre polynomial in x.
        poly_deg_y : int
            Degree of Legendre polynomial in y.
        max_x : int
            Value used to normalize x values.
        '''
        self.max_x = max_x
        # Normalization lambda function
        self.normalize_x = lambda x: x / max_x * 2 - 1
        # x values normalized to [-1, 1]
        self.x_norm = x_norm
        # Physical echelle orders
        self.orders = orders
        # Lambda for normalization of orders
        self.normalize_order = lambda o: (o - np.min(orders)) / (np.max(orders) - np.min(orders)) * 2 - 1
        self.norm_orders = self.normalize_order(orders)
        # Weights
        self.weights = weights
        # Index of lines to include for fitting
        self.index_include = np.array(self.weights>0, dtype=bool)
        # Wavelengths of lines
        self.wavelengths = wavelengths
        # Polynomial degree in x direction
        self.poly_deg_x = poly_deg_x
        # Polynomial degree in y direction
        self.poly_deg_y = poly_deg_y
        # Mean residuals per order
        self.order_means = None

        # 2D Legendre polynomial model
        self.solution = None
        self.residuals = np.zeros_like(self.x_norm)
        self.make_solution()

        # Logger
        self.log = logutils.get_logger(__name__)

    def legendre_fit(self):
        """
        Calculates a 2D Legendre polynomial wavelength fit to normalized x
        and order values.

        Parameters
        ----------
        x_norm : numpy array
            Normalized x values.
        order_norm : numpy array
            Normalized order values.
        weights : numpy array
            Statistical weights for each pixel.
        orders : numpy array
            True order numbers
        wavelengths : numpy array
            Wavelength vectors.
        poly_deg_x : int
            Polynomial degree in x direction
        poly_deg_y : int
            Polynomial degree in y direction

        Returns
        -------
        model : astropy.modeling.Model
            2D wavelength solution polynomial
        """
        log = self.log
        x_norm = self.x_norm
        order_norm = self.norm_orders
        weights = self.weights
        wavelengths = self.wavelengths
        orders = self.orders
        poly_deg_x = self.poly_deg_x
        poly_deg_y = self.poly_deg_y

        init_model = models.Legendre2D(degree=(poly_deg_x, poly_deg_y))
        fit_model = fitting.LinearLSQFitter()
        index_include = np.logical_and(np.array(weights, dtype=bool), np.logical_not(np.isnan(x_norm)))

        t_init = time.time()
        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            model = fit_model(init_model, x_norm[index_include], order_norm[index_include],
                    wavelengths[index_include] * orders[index_include],
                    weights = weights[index_include])
        t_final = time.time()

        log.debug('Time to fit 2D Legendre polynomial: {}'.format(t_final - t_init))
        self.solution = model

    def make_solution(self, weighted=False):
        """
        Make wavelength solution.

        Parameters
        ----------
        weighted : bool
            If true, use weights in the fit.

        Returns
        -------
        None
        """
        if not weighted:
            weights = self.weights[self.index_include] * 0 + 1
        else:
            weights = self.weights[self.index_include]
        self.legendre_fit()

        # Calculate residuals
        self.residuals[self.index_include] = self.solution(self.x_norm[self.index_include],
                                            self.norm_orders[self.index_include]) / self.orders[self.index_include] - \
                                            self.wavelengths[self.index_include]

    def get_wavelength(self, x, order):
        """
        Get wavelength for the specified x value and order.

        Parameters
        ----------
        x : float
            X value.
        order : int
            Order number.

        Returns
        -------
        wavelength : float
            Wavelength for the specified x value and order.
        """
        if isinstance(order, int) or isinstance(order, np.int64):
            order = np.repeat(order, len(x))
        return self.solution(self.normalize_x(x), self.normalize_order(order)) / order

    @deprecated(since="DRAGONS integration", message=_DEPRECATION_MSG)
    def sigma_clip_and_refit(self, threshold=4.0, N=8, weighted=False):
        """
        Sigma clip and refit wavelength solution.

        Parameters
        ----------
        threshold : float
            Sigma clipping threshold.
        N : int
            Number of iterations.
        weighted : bool
            If true, use weights in the fit.

        Returns
        -------
        None
        """
        log = self.log
        residuals = self.residuals[self.index_include]
        wavelengths = self.wavelengths[self.index_include]
        values = residuals/wavelengths * 3e8
        residuals_std = np.std(~np.isnan(values))
        log.info("Start Sigma clipping.  RMS: {:.2f} m/s over {} lines".format(residuals_std, len(values)))
        for i in range(N):
            residuals = self.residuals[self.index_include] # update residuals
            wavelengths = self.wavelengths[self.index_include] # update wavelengths
            values = residuals/wavelengths * 3e8
            residuals_std = np.std(~np.isnan(values))
            if np.count_nonzero(np.abs(values) > threshold * residuals_std) == 0:
                log.info("No more points to clip.  RMS: {:.2f} m/s over {} lines".format(residuals_std, len(values)))
                break
            self.index_include = np.logical_and(self.index_include, np.abs(values) < threshold * residuals_std)
            log.info("Iteration {}: {} points clipped.  RMS: {:.2f} m/s over {} lines".format(i, len(values) - np.count_nonzero(self.index_include), residuals_std, len(values)))
            values = values[self.index_include]
            residuals_std = np.std(~np.isnan(values))
            log.info("RMS after clipping: {:.2f} m/s over {} lines".format(residuals_std, len(values)))
            self.make_solution(weighted=weighted)

    @deprecated(since="DRAGONS integration", message=_DEPRECATION_MSG)
    def calculate_order_means(self):
        '''
        Determines mean of the residuals per order in m/s.

        Args:
            None
        Returns:
            None.  Order means are stored in self.order_means
        '''
         # determines mean of the residuals per order in m/s
        order_means = {}
        unique_orders = np.unique(self.orders)
        wavelengths = np.asarray(self.wavelengths[self.index_include])

        for o in unique_orders:
            residuals = (self.residuals[self.index_include] / wavelengths * 3e8)[
                np.where(self.orders[self.index_include] == o)]
            weights   = (self.weights[self.index_include])[
                np.where(self.orders[self.index_include] == o)]
            order_means[o] = np.average(residuals, weights=weights)
        self.order_means = order_means

