'''
This parameter file contains the parameters related to the primitives
defined in the primitive_maroonx_spectrum.py file.
'''
from gempy.library import config
from geminidr.core import parameters_stack, parameters_standardize
from geminidr.core import parameters_ccd
from geminidr.core import parameters_calibdb  # import core pkgs as needed.

class getPeaksAndPolynomialsConfig(config.Config):
    '''
    This parameter set controls the getPeaksAndPolynomials primitive for MAROON-X.
    '''
    suffix = config.Field("Filename suffix", str, "_wavecal")
    degree_sigma = config.Field("Degree of the sigma polynomial", int, 4)
    degree_width = config.Field("Degree of the width polynomial", int, 2)
    use_sigma_lr = config.Field("Use different polynomials for the left and right sides of the wings", bool, True)
    #plot_path = config.Field("Path to save plots", str, None)
    multithreading = config.Field("Use multithreading", bool, True)
    iterations = config.Field("Number of iterations", int, 3)
