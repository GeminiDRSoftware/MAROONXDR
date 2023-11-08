'''
This module contains the MetaParameter and Parameter classes. The MetaParameter class
is a namedtuple that holds the static meta parameters of the fit and information
derived from those to  handle the parameter vector. The Parameter class is a wrapper
around the parameter vector that allows to easily access the different parts of the
parameter vector and to update it.  Parameter objects also contain the MetaParameter
tuple as a member variable.
'''

from collections import namedtuple
from gempy.utils import logutils
import numpy as np

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
    namedtuple holding the static meta parameters of the fit and information
    derived from those to  handle the parameter vector

    Note:  Immutable
    """

    #Default values
    OFFSET = 0
    SIGMA_LEFT = 1
    SIGMA_RIGHT = 2
    WIDTH = 3
    CENTERS = 4
    AMPLITUDES = 5

    def __new__(cls, sigma, width, number_of_peaks=0, use_sigma_lr=False):
        """
        Create a new MetaParameter object
        Parameters
        ----------
        sigma: int
            number of sigma coefficients
        width: int
            number of width coefficients
        number_of_peaks: int
            number of etalon peaks
        use_sigma_lr: bool
            use different sigmas for left and right side of the wings
        Returns
        -------
        MetaParameter object
        """

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
        Return as numpy array
        Params
        -------
        self

        Returns
        -------
        A numpy array made from the Meta Parameters
        """
        return np.array(self)

    def change_peaks(self, number_of_peaks):
        """
        Return a copy with changed number of peaks

        Params
        -------
        self
        number_of_peaks: int
        A new number of peaks for the meta parameters

        Returns
        -------
        A copy of the Meta Parameters with the new number of peaks
        """
        return MetaParameter(self.sigma, self.width, number_of_peaks, self.use_sigma_lr)

    @property
    def split_at(self):
        """
        Returns the split_at array
        """
        return np.cumsum(
            [
                0,
                self.offset,
                self.sigma_left + 1,
                self.sigma_right + 1,
                self.width + 1,
                self.number_of_peaks,
                self.number_of_peaks,
            ]
        )

    def num_polynomial_parameters(self):
        """
        Returns the number of polynomial parameters (incl. offset)
        """
        return self.split_at[-1]

    def __getnewargs__(self):
        """
        Getter function for the argments of the new meta parameters
        """
        return self.sigma, self.width, self.number_of_peaks, self.use_sigma_lr

    @property
    def offset(self):
        """
        Return the index of the offset
        """
        return self.indices[self.OFFSET]

    @property
    def centers(self):
        """
        Return the indices of the centers of the etalon peaks
        """
        return self.indices[self.CENTERS]

    @property
    def amplitudes(self):
        """
        Return the indices of the amplitudes of the etalon peaks
        """
        return self.indices[self.AMPLITUDES]

    @property
    def polynomials(self):
        """
        Return the indices of the polynomials
        """
        return self.indices[self.SIGMA_LEFT: -2]

class Parameter(object):
    '''
    Wrapper around the parameter vector that allows to easily access the different
    parts of the parameter vector and to update it.  Parameter objects also contain
    the MetaParameter tuple as a member variable.
    '''
    def __init__(self, offset = 0,
                 p_sigma_left = None,
                 p_sigma_right = None,
                 p_width = None,
                 peaks = None,
                 amplitudes = None,
                 meta_parameters = None):
        """
        Concatenate a parameter vector into a single numpy array to pass it to
        scipy.least_squares

        Parameters
        ----------
        offset: int
            offset
        p_sigma_left: list of ints
            sigmas to use for the left side
        p_sigma_right: list of ints
            sigmas to use for the right side
        p_width: list of ints
            coefficients of width polynomial
        peaks: list of ints
            center of etalon peaks
        amplitudes: list of ints
            amplitudes of peaks
        meta_parameters: MetaParameters
            Fit meta parameters

        Returns
        -------
        res : ndarray
            concatenated fit parameters
        """

        self.log = logutils.get_logger(__name__)

        if not meta_parameters.use_sigma_lr:
            assert np.all(p_sigma_left == p_sigma_right)
            p_sigma_right = np.array([])

        res = np.concatenate([[offset], p_sigma_left, p_sigma_right, p_width, peaks, amplitudes])
        assert len(res) == meta_parameters.total, f"{len(res)} != {meta_parameters.total}"
        self.parameters = res
        self.meta_parameters = meta_parameters

    @property
    def offset(self):
        """
        Get the offset from the parameters.
        Overloaded function that can
        """
        return self.parameters[0] # offset is always the first parameter

    @property
    def centers(self):
        """
        Get the centers of the etalon peaks from the parameters
        """
        return self.parameters[self.meta_parameters.centers]

    @property
    def amplitudes(self):
        """
        Get the amplitudes of the etalon peaks from the parameters
        """
        return self.parameters[self.meta_parameters.amplitudes]

    def eval_polynomials(self, x):
        """
        Evaluate the polynomials for the given parameters
        Parameters
        ----------
        self : Parameter object
        x: ndarray
            x values to evaluate the polynomials at
        Returns
        -------
        values: ndarray
            values of the polynomials at the given x values
        """
        polynomials = self.meta_parameters.polynomials
        values = np.zeros(shape=(len(polynomials), len(x)))
        for ii, idx in enumerate(polynomials):
            values[ii] = np.poly1d(self.parameters[idx])(x)
        return values

    def eval_polynomials_at_centers(self):
        """
        Evaluate the polynomials at the centers of the fitted etalon peaks
        Parameters
        ----------
        parameters: ndarray
            parameters of the polynomials

        Returns
        -------
        ndarray
            values of the polynomials at the centers of the
            fitted etalon peaks
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
        Split parameter array into subarrays, reverses concat_parameters

        Parameters
        ----------
        parameters: ndarray
            parameter vector
        meta_parameters: MetaParameters object
            fit meta parameters
        Returns
        ---------
        values: tuple
        (offset, sigma coefficients,
                width coefficient, peak center locations, amplitudes)
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
        Returns a copy of the parameter vector with updated parameters
        Parameters
        ----------
        parameters: ndarray
            parameter vector
        meta_parameters: MetaParameters object
            fit meta parameters
        offset: int
            offset
        p_sigma_l: list of ints
            sigmas to use for the left side
        p_sigma_r: list of ints
            sigmas to use for the right side
        p_width: list of ints
            coefficients of width polynomial
        centers: list of ints
            center of etalon peaks
        amplitudes: list of ints
            amplitudes of peaks
        Returns
        -------
        p : ndarray
            updated parameter vector
        """
        p = self.parameters.copy()
        values = [offset, p_sigma_l, p_sigma_r, p_width, centers, amplitudes]
        for idx, v in zip(self.meta_parameters.indices, values):
            if v is not None:
                p[idx] = v
        if parameters is not None:
            self.parameters = parameters.copy()
