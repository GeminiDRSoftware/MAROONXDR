"""
This parameter file contains the parameters related to the primitives
defined in the primitive_maroonx_spectrum.py file.
"""

from gempy.library import config


class staticWavelengthSolutionConfig(config.Config):
    """
    This parameter set controls the staticWavelengthSolution primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')
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

    suffix = config.Field('Filename suffix', str, '')
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
    suffix = config.Field('Filename suffix', str, '_etalonwls')
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
    report = config.Field(
        'Write PDF report with diagnostic plots',
        bool,
        True,
    )


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
    report = config.Field(
        'Write PDF report with diagnostic plots',
        bool,
        True,
    )

class combineFibersConfig(config.Config):
    """
    This parameter set controls the combineFibers primitive for MAROON-X.
    """

    suffix = config.Field('Filename suffix', str, '')
    combine_fibers = config.ListField(
        'List of fibers to combine.',
        int,
        None,
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
    report = config.Field('Generate PDF diagnostic report', bool, True)


class barycentricCorrectionConfig(config.Config):
    """
    This parameter set controls the barycentricCorrection primitive for MAROON-X.
    """
    suffix = config.Field('Filename suffix', str, '_reduced')
    target_name = config.Field(
        'Target name to downselect files',
        str,
        None,
        optional=True,
    )
    simbad_target_name = config.Field(
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
    report = config.Field(
        "Generate diagnostic PDF of exposure meter time series.",
        bool,
        True,
    )


class separateArmStreamsConfig(config.Config):
    """
    Configuration for separateArmStreams primitive.
    """
    pass


class bundleArmStreamsConfig(config.Config):
    """
    Configuration for bundleArmStreams primitive.
    """
    suffix = config.Field('Filename suffix', str, '_rebundled')


class displaySpectraConfig(config.Config):
    """
    Configuration for displaySpectra primitive.

    This primitive launches an interactive Bokeh viewer to display extracted
    spectra in a web browser for quality assessment.
    """
    fibers = config.ListField(
        'List of fibers to display (e.g., [2, 3, 4] for science fibers, '
        '[6] for combined fiber, [5] for calibration fiber).',
        int,
        None,
        optional=True,
        single=False,
    )
    show_wavelength = config.Field(
        'Display spectra vs wavelength (nm) if wavelength solution available. '
        'If False or no wavelength solution exists, display vs pixel number.',
        bool,
        False,
    )