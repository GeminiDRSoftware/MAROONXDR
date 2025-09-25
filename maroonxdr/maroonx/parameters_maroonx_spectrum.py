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

    suffix = config.Field('Filename suffix', str, '_peaks')
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
    orders = config.ListField(
        'List of orders to process.', 
        int, 
        None, 
        optional=True, 
        single=True,
    )
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
    fibers = config.ListField(
        'List of fibers to process.',
        int,
        None,
        optional=True,
        single=True,
    )
    symmetric_linefits = config.Field('Symmetric line fits', bool, False)
    n_knots = config.Field(
        'Number of knots for the cubic spline fit',
        int,
        30,
    )
    thar = config.Field(
        'Whether to apply ThAr wavelength solution to Etalon frames',
        bool,
        False,
    )
    ref_file = config.Field('Reference file', str, None, optional=True)
    ref_fiber = config.Field('Reference fiber', int, None, optional=True)
    #plot_path = config.Field('Path to save plots', str, '')


class applyWavelengthSolutionConfig(config.Config):
    """
    This parameter set controls the applyWavelengthSolution primitive for MAROON-X.
    """
    suffix = config.Field('Filename suffix', str, '_wls')
    fibers = config.ListField(
        'List of fibers to process.',
        int,
        None,
        optional=True,
        single=True,
    )
    symmetric_linefits = config.Field('Symmetric line fits', bool, False)
    n_knots = config.Field(
        'Number of knots for the cubic spline fit',
        int,
        30,
    )
    thar = config.Field(
        'Whether to apply ThAr wavelength solution to Etalon frames',
        bool,
        False,
    )
    ref_fiber = config.Field('Reference fiber', int, None, optional=True)
    #plot_path = config.Field('Path to save plots', str, '')

class combineFibersConfig(config.Config):
    """
    This parameter set controls the combineFibers primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '_combined')
    combine_fibers = config.ListField(
        'List of fibers to combine.',
        int,
        [2, 3, 4],
        optional=True,
        single=False,
    )
    symmetric_linefits = config.Field('Symmetric line fits', bool, False)
    kappa_sigma = config.Field(
        'Sigma clipping threshold for outlier rejection',
        float,
        5.,
        optional=True,
    )
    max_clips = config.Field(
        'Maximum pixels to clip per order before increasing kappa_sigma',
        int,
        5000,
        optional=True,
    )

class barycentricCorrectionConfig(config.Config):
    """
    This parameter set controls the barycentricCorrection primitive for MAROON-X.
    """
    suffix = config.Field('Filename suffix', str, '_reduced')
    target_name = config.ListField(
        'Target name to downselect files',
        str,
        None,
        optional=True,
    )
    simbad_target_name = config.ListField(
        'SIMBAD resolvable target name',
        str,
        None,
        optional=True,
    )
    use_coords = config.Field(
        'Use telescope pointing coordinates instead of target name',
        bool,
        False,
    )
    zp_pc = config.Field(
        'Zeropoint for counts_pc channel. Determined from data if not provided.',
        float,
        0.,
        optional=True,
    )
    zp_frd = config.Field(
        'Zeropoint for counts_frd channel. Determined from data if not provided.',
        float,
        0.,
        optional=True,
    )
    start_time = config.ChoiceField(
        "Time to consider to compute exposure start",
        str,
        allowed={"mjd_start": "Telescope MJD written at start of exposure",
                 "mjd_end": "Telescope MJD written at end of readout",
                 "filename": "UTC from filename"},
        default="filename")