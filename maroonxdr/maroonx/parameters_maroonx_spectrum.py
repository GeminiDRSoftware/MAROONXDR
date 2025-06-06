"""
This parameter file contains the parameters related to the primitives
defined in the primitive_maroonx_spectrum.py file.
"""

from gempy.library import config


class staticWavelengthSolutionConfig(config.Config):
    """
    This parameter set controls the staticWavelengthSolution primitive for MAROON-X.
    """

    fibers = config.ListField(
        'List of fibers to process.',
        int,
        None,
        optional=True,
        single=True,
    )


class getPeaksAndPolynomialsConfig(config.Config):
    """
    This parameter set controls the getPeaksAndPolynomials primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '_wavecal')
    guess_file = config.Field(
        'Name of file containing initial guess spectrum.',
        str,
        None,
        optional=True
    )
    fibers = config.ListField(
        'List of fibers to process.',
        int,
        None,
        optional=True,
        single=True,
    )
    orders = config.ListField('Orders to fit.', int, None, optional=True)
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

    #plot_path = config.Field('Path to save plots', str, '')
    ref_file = config.Field('Reference file', str, None, optional=True)
    ref_fiber = config.Field('Reference fiber', int, 5)
    symmetric_linefits = config.Field('Symmetric line fits', bool, False)
