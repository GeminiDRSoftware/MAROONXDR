from .echellespectrum import EchelleSpectrum

'''
This class is for ThAr Spectra. Debatable whether we will want to keep this file, as static wavelength solutions
are calculated extremely rarely and do not need to be part of Gemini's Data Reduction Pipeline.
'''
class ThArSpectrum(EchelleSpectrum):
    """
    Describes a ThAr spectrum.  Inherits from EchelleSpectrum.
    """
    def __init__(self, **kwargs):
        """
        Initializes the ThArSpectrum object.
        """
        pass

    def load_manual_guesses(self, x, wavelengths, orders, auto_fit=True, auto_cc=True, fit_range=4, debug_plots = False):
        """
        Load manual guesses pixel->wavelength per order.

        Parameters
        ----------
        debug : bool
            Debug flag.  If True, plots will be generated.
        fit_range : int
            Number of pixels on either side of the peak to use for fitting.
        auto_fit : bool
            If true, manual guesses are re-fitted.
            TODO:  Use Peak Modeller for this as opposed to Gaussian Model.
        x: Iterable
            Guessed x-position of a line. Same length as orders.
        wavelengths: Iterable
            Catalogue vaccuum wavelengths of the lines. Same length as orders.
        orders: Iterable
            Vector of physical echelle orders. Same length as x and wavelengths.
        auto_cc: bool
            If true, manual guesses are cross-correlated with the data to calculate the shift between
            the guesses and the data.
        """
        pass

    def plot_wavelength_guess(self):
        """
        Plot manual wavelength guesses on spectrume

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        pass

    def fit_initial_lines(self, data, thar_pos_guess, fiber, order, fit_range=5, debug=0):
        """
        Fit Thorium Argon emission lines from initial (manual guess.
        It uses the spectrum's peak model to fit the peaks.

        Parameters
        ----------
        data : 1D array
            Data to fit.
        thar_pos_guess : 1D array
            Initial guesses for the peak positions.
        fiber : int
            Fiber number.
        order : int
            Order number.
        fit_range : int
            Number of pixels on either side of the peak to use for fitting.
        debug : int
            Debug level.
        Returns
        -------
        Fitted ThAr peaks
        """
        pass

    def make_initial_wavelength_solutions(self, poly_deg_x=5, poly_deg_y=5, debug=0):
        """
        Make initial wavelength solutions for each order.

        Parameters
        ----------
        poly_deg_x : int
            Polynomial degree in x.
        poly_deg_y : int
            Polynomial degree in y.
        debug : int
            Debug level.

        Returns
        -------
        None
        """
        pass
