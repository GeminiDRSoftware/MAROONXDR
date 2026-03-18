import os
import warnings
from pathlib import Path
from urllib.error import HTTPError

import numpy as np
import pytest
from astrodata.testing import download_from_archive

# =========================================================
# DRAGONS TEST CONFIGURATION
# =========================================================


def get_dragons_test_path():
    """
    Get the DRAGONS_TEST environment variable path.

    Returns
    -------
    Path or None
        Path to DRAGONS test data directory, or None if not set
    """
    p = os.environ.get('DRAGONS_TEST')
    return Path(p) if p else None


def get_maroonx_legacy_test_path():
    """
    Get the legacy environment variable path to all
    maroonx reduced files using the legacy pipeline.

    Returns
    -------
    Path or None
        Path to legacy test data directory, or None if not set
    """
    p = os.environ.get('MAROONX_LEGACY_TEST')
    return Path(p) if p else None


# =========================================================
# PATH FIXTURES
# =========================================================

# Dragons related fixtures

@pytest.fixture(scope='session')
def dragons_test_root():
    """
    Fixture providing the root DRAGONS test data directory.

    Returns
    -------
    Path
        Root path to DRAGONS test data
    """
    root = get_dragons_test_path()
    if root is None:
        pytest.skip('DRAGONS_TEST environment variable not set')
    if not root.exists():
        pytest.skip(f'DRAGONS test root does not exist: {root}')
    return root


@pytest.fixture(scope='session')
def preprocessed_files_path(dragons_test_root):
    """Path to preprocessed (debundled) FITS files."""
    path = dragons_test_root / 'preprocessed_files'
    if not path.exists():
        pytest.skip(f'Preprocessed files directory does not exist: {path}')
    return path


@pytest.fixture(scope='function')
def processed_dark_path(dragons_test_root):
    """Path to processed dark calibrations."""
    path = dragons_test_root / 'preprocessed_files' / 'calibrations' / 'processed_dark'
    if not path.exists():
        pytest.skip(f'Processed dark directory does not exist: {path}')
    return path


@pytest.fixture(scope='function')
def processed_dark_coeff_path(dragons_test_root):
    """Path to processed dark coefficient files."""
    path = dragons_test_root / 'preprocessed_files' / 'calibrations' / 'processed_dark_coeff'
    if not path.exists():
        pytest.skip(f'Processed dark coeff directory does not exist: {path}')
    return path


# Legacy pipeline fixtures

@pytest.fixture(scope='session')
def legacy_test_root():
    """
    Fixture providing the root legacy test data directory.

    Returns
    -------
    Path
        Root path to legacy test data
    """
    root = get_maroonx_legacy_test_path()
    if root is None:
        pytest.skip('MAROONX_LEGACY_TEST environment variable not set')
    if not root.exists():
        pytest.skip(f'Legacy root does not exist: {root}')
    return root


@pytest.fixture(scope='function')
def path_to_legacy_darks(legacy_test_root):
    """
    Fixture providing path to test legacy data directory.
    """
    path = (
        legacy_test_root
        / 'MaroonX_spectra_reduced'
        / 'Maroonx_masterframes'
        / '202411xx'
        / 'darks'
    )
    if not path.exists():
        pytest.skip(f'Legacy data directory does not exist: {path}')
    return path


@pytest.fixture(scope='function')
def path_to_legacy_flats(legacy_test_root):
    """
    Fixture providing path to test legacy data directory.
    """
    path = (
        legacy_test_root
        / 'MaroonX_spectra_reduced'
        / 'Maroonx_masterframes'
        / '202411xx'
        / 'flats'
    )
    if not path.exists():
        pytest.skip(f'Legacy data directory does not exist: {path}')
    return path


@pytest.fixture(scope='function')
def path_to_legacy_wavecal(legacy_test_root):
    """
    Fixture providing path to test legacy data directory.
    """
    path = legacy_test_root / 'MaroonX_spectra_reduced' / '20241124'
    if not path.exists():
        pytest.skip(f'Legacy data directory does not exist: {path}')
    return path


@pytest.fixture(scope='function')
def path_to_legacy_science(legacy_test_root):
    """
    Fixture providing path to test legacy data directory.
    """
    path = legacy_test_root / 'MaroonX_spectra_reduced' / '20241124'
    if not path.exists():
        pytest.skip(f'Legacy data directory does not exist: {path}')
    return path


@pytest.fixture(scope='function')
def path_to_legacy_reduced(legacy_test_root):
    """
    Fixture providing path to legacy reduced science data.

    Contains both intermediate .hdf files and final .hd5 files
    produced by the legacy MaroonX pipeline.
    """
    path = legacy_test_root / 'MaroonX_spectra_reduced' / '20241124'
    if not path.exists():
        pytest.skip(f'Legacy reduced data directory does not exist: {path}')
    return path


@pytest.fixture(scope='session')
def path_to_legacy_bkg(legacy_test_root):
    """Path to legacy background/intermediate npy arrays."""
    path = legacy_test_root.parent / 'legacy_bkg_arrays'
    if not path.exists():
        pytest.skip(f'Legacy bkg arrays directory does not exist: {path}')
    return path


# =========================================================
# FIXTURES
# =========================================================


@pytest.fixture(params=['BLUE', 'RED'], scope='function')
def arm(request):
    """
    Fixture providing color tag name for the arm.
    """
    return request.param


# =========================================================
# UTILITY FUNCTIONS
# =========================================================


def assert_allclose_with_max_fails(
    x, y, rtol, atol, max_fails=0, warn_on_acceptable_fails=True
):
    """
    Assert arrays are close, allowing up to max_fails elements to fail.

    This is useful for regression testing where a small number of pixel-level
    differences are acceptable due to rounding errors or algorithm variations.

    When max_fails > 0 and failures are within acceptable range (0 < n_fails <= max_fails),
    the test will trigger a pytest.xfail() to mark it as an expected failure rather than
    a clean pass. This provides better visibility into which tests have acceptable
    deviations versus perfect matches.

    Parameters
    ----------
    x : array_like
        First array to compare
    y : array_like
        Second array to compare
    rtol : float
        Relative tolerance parameter. No default, must be specified.
    atol : float
        Absolute tolerance parameter. No default, must be specified.
    max_fails : int, optional
        Maximum number of elements allowed to fail (default: 0)
    warn_on_acceptable_fails : bool, optional
        If True, trigger pytest.xfail() when 0 < n_fails <= max_fails (default: True)
        If False, silently pass when failures are within acceptable range

    Raises
    ------
    AssertionError
        If more than max_fails elements fail the tolerance test
    """
    x = np.asarray(x)
    y = np.asarray(y)

    # Check which elements fail the tolerance test
    not_close = ~np.isclose(x, y, rtol=rtol, atol=atol)
    n_fails = np.sum(not_close)

    if n_fails > max_fails:
        fail_indices = np.where(not_close)[0]
        raise AssertionError(
            f'{n_fails} elements failed (max allowed: {max_fails})\n'
            f'First few failing indices: {fail_indices[:10]}'
        )
    if n_fails > 0 and warn_on_acceptable_fails:
        # Mark as expected failure
        fail_indices = np.where(not_close)[0]
        pytest.xfail(
            f'{n_fails} elements differ (within acceptable limit of {max_fails})\n'
        )


# =========================================================
# CONFIGURATION HOOKS
# =========================================================


def pytest_addoption(parser):
    """
    Add custom command-line options for MaroonX tests.

    Parameters
    ----------
    parser : pytest.Parser
        Pytest argument parser
    """
    parser.addoption(
        '--preprocess-bundles',
        action='store_true',
        default=False,
        help='Run bundle preprocessing before tests (splits bundles)',
    )
    parser.addoption(
        '--create-inputs',
        action='store_true',
        default=False,
        help='Create test input files in $DRAGONS_TEST directory structure',
    )


def pytest_configure(config):
    """
    Configure pytest with custom markers and settings.

    Parameters
    ----------
    config : pytest.Config
        Pytest configuration object
    """
    # Register custom markers
    config.addinivalue_line(
        'markers', 'maroonx: marks tests as MAROON-X instrument tests'
    )
    config.addinivalue_line(
        'markers',
        'regression: marks tests as regression tests '
        '(deselect with \'-m "not regression"\')',
    )
    config.addinivalue_line(
        'markers',
        'legacy_regression: marks tests as legacy pipeline regression tests '
        '(deselect with \'-m "not legacy_regression"\')',
    )
    config.addinivalue_line(
        'markers', 'slow: mark test as slow tests (deselect with \'-m "not slow"\')'
    )
    config.addinivalue_line(
        'markers',
        'preprocessed_data: marks tests that require preprocessed data in inputs/ '
        '(deselect with \'-m "not preprocessed_data"\')',
    )


def pytest_collection_modifyitems(config, items):
    """
    Modify test collection to add automatic markers based on test location and names.

    Parameters
    ----------
    config : pytest.Config
        Pytest configuration object
    items : list
        List of collected test items
    """
    for item in items:
        # Auto-apply maroonx marker to all tests in maroonx test suite
        if 'maroonx/tests' in str(item.path):
            item.add_marker(pytest.mark.maroonx)

        # Auto-apply regression marker to tests in regression directories
        if 'regression' in str(item.path):
            item.add_marker(pytest.mark.regression)

        # Auto-apply legacy_regression marker to tests in legacy_regression directory
        if 'legacy_regression' in str(item.path):
            item.add_marker(pytest.mark.legacy_regression)


# =========================================================
# MANIFEST AND ARCHIVE DOWNLOAD
# =========================================================

# Manifest of MaroonX raw files available in the Gemini Archive
# These files were used for MaroonX DRAGONS development and testing
MAROONX_TEST_MANIFEST = {
    'DARK': [
        'N20241115M3421.fits',
        'N20241115M3433.fits',
        'N20241115M3444.fits',
        'N20241115M3455.fits',
        'N20241115M3466.fits',
        'N20241115M3477.fits',
        'N20241115M3486.fits',
        'N20241115M3494.fits',
        'N20241115M3502.fits',
        'N20241115M3510.fits',
        'N20241115M3519.fits',
        'N20241115M3539.fits',
        'N20241115M3559.fits',
        'N20241115M3579.fits',
        'N20241115M3600.fits',
        'N20241115M3620.fits',
        'N20241115M3655.fits',
        'N20241115M3690.fits',
        'N20241115M3726.fits',
        'N20241115M3761.fits',
        'N20241115M3796.fits',
        'N20241115M3846.fits',
        'N20241115M3897.fits',
        'N20241115M3947.fits',
        'N20241115M3997.fits',
        'N20241115M4047.fits',
        'N20241115M4112.fits',
        'N20241115M4178.fits',
        'N20241115M4243.fits',
        'N20241115M4308.fits',
        'N20241116M0054.fits',
        'N20241116M0149.fits',
        'N20241116M0244.fits',
        'N20241116M0339.fits',
        'N20241116M0434.fits',
        'N20241124M0655.fits',
        'N20241124M0659.fits',
        'N20241124M0663.fits',
        'N20241124M0668.fits',
        'N20241124M0672.fits',
        'N20241124M3038.fits',
        'N20241124M3043.fits',
        'N20241124M3047.fits',
        'N20241124M3051.fits',
        'N20241124M3055.fits',
        'N20241125M0774.fits',
    ],
    'FLAT': [
        'N20241114M3271.fits',
        'N20241114M3290.fits',
        'N20241114M3295.fits',
        'N20241114M3300.fits',
        'N20241114M3305.fits',
        'N20241114M3310.fits',
        'N20241114M3442.fits',
        'N20241114M3450.fits',
        'N20241114M3456.fits',
        'N20241114M3461.fits',
        'N20241114M3466.fits',
        'N20241114M3471.fits',
    ],
    'WAVECAL': [
        'N20240814M0349.fits',
        'N20241124M0542.fits',
        'N20241124M0547.fits',
        'N20241124M0554.fits',
        'N20241124M0559.fits',
        'N20241124M0617.fits',
        'N20241124M0622.fits',
        'N20241124M0639.fits',
        'N20241124M2945.fits',
        'N20241124M2951.fits',
        'N20241124M2957.fits',
        'N20241124M2962.fits',
        'N20241124M3018.fits',
        'N20241124M3023.fits',
        'N20241124M3032.fits',
    ],
    'SCIENCE': [
        # "N20241124M1116.fits",
    ],
}


@pytest.fixture(scope='session')
def download_mx_file(dragons_test_root):
    """
    Session fixture that provides a function to download MaroonX files from the Gemini Archive.

    If a file returns HTTP 403 (Access Forbidden), the test will be skipped.
    """

    def _download(filename, sub_path='science_dir'):
        try:
            return download_from_archive(
                filename, sub_path=sub_path, env_var='DRAGONS_TEST'
            )
        except HTTPError as e:
            # dont fail if one file is not accessible
            warnings.warn(f'{filename}: {e}')
            return None

    return _download
