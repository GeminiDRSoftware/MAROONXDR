"""
Parameter-vector handling for the etalon spectrum fit.

``MetaParameter`` is an immutable namedtuple holding the static meta
parameters of the fit and the derived slices that select the parts of
the flat parameter vector. ``Parameter`` wraps the parameter vector itself,
giving named access to its parts and in-place updates, and carries the
``MetaParameter`` tuple as a member variable.
"""

from collections import namedtuple
import numpy as np

from . import get_logger


MetaParameterBase = namedtuple(
        "MetaParameterBase",
        [
            "offset",
            "sigma",
            "width",
            "number_of_peaks",
            "total",
            "indices",
            "use_sigma_lr",
        ],
)


class MetaParameter(MetaParameterBase):
    """
    Immutable namedtuple holding the static meta parameters of the fit.

    Also holds the information derived from those to select the parts
    of the flat parameter vector.

    Parameters
    ----------
    sigma : int
        Degree of the sigma polynomial(s).

    width : int
        Degree of the width polynomial.

    number_of_peaks : int
        Number of etalon peaks. Default is 0.

    use_sigma_lr : bool
        Use different sigmas for the left and right flanks of the peaks.
        Default is False.

    Attributes
    ----------
    total : int
        Total length of the parameter vector.

    indices : list of slice
        Slices that select the parts of the parameter vector, in order:
        offset, left sigma coefficients, right sigma coefficients, width
        coefficients, peak centers, amplitudes. When ``use_sigma_lr`` is
        False the right sigma slice equals the left one.
    """

    #Default values
    OFFSET = 0
    SIGMA_LEFT = 1
    SIGMA_RIGHT = 2
    WIDTH = 3
    CENTERS = 4
    AMPLITUDES = 5

    def __new__(cls, sigma, width, number_of_peaks=0, use_sigma_lr=False):
        offset = 1
        sigma_left = sigma
        sigma_right = sigma if use_sigma_lr else -1
        # Generate a list of indices to split the parameter vector
        split_at = np.cumsum(
            [
                0,
                offset,
                sigma_left + 1,
                sigma_right + 1,
                width + 1,
                number_of_peaks,
                number_of_peaks,
            ]
        )

        # Create a list of slices from the split_at array to actually allow us to split the vector
        indices = [slice(l, r) for l, r in zip(split_at[:-1], split_at[1:])]

        if not use_sigma_lr:
            # Use the same values for both sigmas- we can do this by making those 2 slices the same
            indices[2] = indices[1]

        total = indices[-1].stop
        self = super().__new__(
            cls, offset, sigma, width, number_of_peaks, total, indices, use_sigma_lr
        )
        return self

    def as_array(self):
        """
        Return the meta parameters as a numpy array.

        Returns
        -------
        ndarray
            Array made from the meta parameter fields.
        """
        return np.array(self)

    def change_peaks(self, number_of_peaks):
        """
        Return a copy with a changed number of peaks.

        Parameters
        ----------
        number_of_peaks : int
            New number of peaks.

        Returns
        -------
        MetaParameter
            Copy of the meta parameters with the new number of peaks.
        """
        return MetaParameter(self.sigma, self.width, number_of_peaks, self.use_sigma_lr)

    def __getnewargs__(self):
        """Return the constructor arguments, used for pickling."""
        return self.sigma, self.width, self.number_of_peaks, self.use_sigma_lr

    @property
    def offset(self):
        """Return the slice that selects the offset from the parameter vector."""
        return self.indices[self.OFFSET]

    @property
    def centers(self):
        """Return the slice that selects the peak centers from the parameter vector."""
        return self.indices[self.CENTERS]

    @property
    def amplitudes(self):
        """Return the slice that selects the peak amplitudes from the parameter vector."""
        return self.indices[self.AMPLITUDES]

    @property
    def polynomials(self):
        """Return the slices that select the polynomial coefficients (left sigma, right sigma, width)."""
        return self.indices[self.SIGMA_LEFT: -2]


class Parameter(object):
    """
    Wrapper around the flat parameter vector of the fit.

    Gives named access to the different parts of the parameter vector
    and allows updating it in place. The constructor concatenates the
    parts into a single numpy array as expected by
    ``scipy.optimize.least_squares``. When
    ``meta_parameters.use_sigma_lr`` is False, ``p_sigma_left`` and
    ``p_sigma_right`` must be equal and only one copy is stored.

    Parameters
    ----------
    offset : float
        Constant offset of the spectrum.

    p_sigma_left : ndarray
        Coefficients of the left sigma polynomial.

    p_sigma_right : ndarray
        Coefficients of the right sigma polynomial.

    p_width : ndarray
        Coefficients of the width polynomial.

    peaks : ndarray
        Centers of the etalon peaks.

    amplitudes : ndarray
        Amplitudes of the peaks.

    meta_parameters : MetaParameter
        Fit meta parameters.

    Attributes
    ----------
    parameters : ndarray
        Concatenated fit parameter vector.

    meta_parameters : MetaParameter
        Fit meta parameters.
    """
    def __init__(self, offset = 0,
                 p_sigma_left = None,
                 p_sigma_right = None,
                 p_width = None,
                 peaks = None,
                 amplitudes = None,
                 meta_parameters = None):
        # self.log = logutils.get_logger(__name__)
        self.log = get_logger()

        if not meta_parameters.use_sigma_lr:
            assert np.all(p_sigma_left == p_sigma_right)
            p_sigma_right = np.array([])

        res = np.concatenate([[offset], p_sigma_left, p_sigma_right, p_width, peaks, amplitudes])
        assert len(res) == meta_parameters.total, f"{len(res)} != {meta_parameters.total}"
        self.parameters = res
        self.meta_parameters = meta_parameters

    @property
    def offset(self):
        """Return the offset from the parameter vector."""
        return self.parameters[0] # offset is always the first parameter

    @property
    def centers(self):
        """Return the centers of the etalon peaks from the parameter vector."""
        return self.parameters[self.meta_parameters.centers]

    @property
    def amplitudes(self):
        """Return the amplitudes of the etalon peaks from the parameter vector."""
        return self.parameters[self.meta_parameters.amplitudes]

    def eval_polynomials(self, x):
        """
        Evaluate the polynomials at the given x values.

        Parameters
        ----------
        x : ndarray
            x values to evaluate the polynomials at.

        Returns
        -------
        ndarray
            Values of the polynomials at the given x values, one row per
            polynomial.
        """
        polynomials = self.meta_parameters.polynomials
        values = np.zeros(shape=(len(polynomials), len(x)))
        for ii, idx in enumerate(polynomials):
            values[ii] = np.poly1d(self.parameters[idx])(x)
        return values

    def eval_polynomials_at_centers(self):
        """
        Evaluate the polynomials at the centers of the fitted etalon peaks.

        Returns
        -------
        ndarray
            Peak centers, amplitudes, and polynomial values at the centers,
            one row each.
        """
        values = [
            self.parameters[self.meta_parameters.centers],
            self.parameters[self.meta_parameters.amplitudes],
        ]
        for idx in self.meta_parameters.polynomials:
            values.append(np.poly1d(self.parameters[idx])(values[0]))
        return np.array(values)

    def split_parameters(self):
        """
        Split the parameter vector into its parts.

        Reverses the concatenation done by the constructor.

        Returns
        -------
        list
            ``[offset, left sigma coefficients, right sigma coefficients,
            width coefficients, peak centers, amplitudes]``, with the
            offset as a scalar.
        """

        values = [self.parameters[idx] for idx in self.meta_parameters.indices]
        centers = values[-2]
        amplitudes = values[-1]
        assert len(centers) == self.meta_parameters.number_of_peaks,\
        f"{centers} != {self.meta_parameters.number_of_peaks}"
        assert len(amplitudes) == self.meta_parameters.number_of_peaks,\
        f"{centers} != {self.meta_parameters.number_of_peaks}"
        values[0] = np.ndarray.item(values[0])
        return values

    def update_parameters(
            self,
            parameters=None,
            offset=None,
            p_sigma_l=None,
            p_sigma_r=None,
            p_width=None,
            centers=None,
            amplitudes=None,
    ):
        """
        Update the parameter vector in place.

        Parts whose corresponding argument is None are left unchanged.
        If ``parameters`` is given, it replaces the whole vector and
        takes precedence.

        Parameters
        ----------
        parameters : ndarray or None
            Complete replacement parameter vector; stored as a copy.

        offset : float or None
            Constant offset of the spectrum.

        p_sigma_l : ndarray or None
            Coefficients of the left sigma polynomial.

        p_sigma_r : ndarray or None
            Coefficients of the right sigma polynomial.

        p_width : ndarray or None
            Coefficients of the width polynomial.

        centers : ndarray or None
            Centers of the etalon peaks.

        amplitudes : ndarray or None
            Amplitudes of the peaks.
        """
        values = [offset, p_sigma_l, p_sigma_r, p_width, centers, amplitudes]
        for idx, v in zip(self.meta_parameters.indices, values):
            if v is not None:
                self.parameters[idx] = v
        if parameters is not None:
            self.parameters = parameters.copy()