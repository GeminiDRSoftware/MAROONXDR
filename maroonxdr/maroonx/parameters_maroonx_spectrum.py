"""
This parameter file contains the parameters related to the primitives
defined in the primitive_maroonx_spectrum.py file.
"""

from gempy.library import config


class getPeaksAndPolynomialsConfig(config.Config):
    """
    This parameter set controls the getPeaksAndPolynomials primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '_wavecal')
    degree_sigma = config.Field('Degree of the sigma polynomial', int, 4)
    degree_width = config.Field('Degree of the width polynomial', int, 2)
    use_sigma_lr = config.Field(
        'Use different polynomials for the left and right sides of the wings',
        bool,
        True,
    )
    # plot_path = config.Field("Path to save plots", str, None)
    multithreading = config.Field('Use multithreading', bool, False)
    iterations = config.Field('Number of iterations', int, 10)


class fitAndApplyEtalonWlsConfig(config.Config):
    """
    This parameter set controls the fitAndApplyEtalonWls primitive for MAROON-X.
    """

    plot_path = config.Field('Path to save plots', str, '')
    ref_file = config.Field('Reference file', str, '')
    ref_fiber = config.Field('Reference fiber', int, 5)
    symmetric_linefits = config.Field('Symmetric line fits', bool, False)
