import os
import warnings
from pathlib import Path
from urllib.error import HTTPError

import numpy as np
import pytest
from astropy.io import fits
from astrodata.testing import download_from_archive

import astrodata
import maroonx_instruments  # noqa: F401 — registers MaroonX AstroData class

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
        / '202507xx'
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
        / '202507xx'
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
    path = legacy_test_root / 'MaroonX_spectra_reduced' / '20250717'
    if not path.exists():
        pytest.skip(f'Legacy data directory does not exist: {path}')
    return path


@pytest.fixture(scope='function')
def path_to_legacy_science(legacy_test_root):
    """
    Fixture providing path to test legacy data directory.
    """
    path = legacy_test_root / 'MaroonX_spectra_reduced' / '20250717'
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
    path = legacy_test_root / 'MaroonX_spectra_reduced' / '20250717'
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


@pytest.fixture(params=['RED', 'BLUE'])
def ad_min(request):
    """Minimal MaroonX AstroData object for synthetic unit tests.

    Uses 4400x4400 arrays matching the real detector size so that lookup-based
    section descriptors (overscan, data, array) work without mocking.  Headers
    are the minimum set needed for tag resolution and descriptor lookup (ARM,
    INSTRUME, FIBERn, EXPTIME).
    """
    arm = request.param
    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('FIBER1', 'Dark')
    phu.header.set('FIBER2', 'Dark')
    phu.header.set('FIBER3', 'Dark')
    phu.header.set('FIBER4', 'Dark')
    phu.header.set('FIBER5', 'Etalon')

    sci = fits.ImageHDU(data=np.ones((4400, 4400), dtype=np.float32), name='SCI')
    sci.header.set('ARM', arm)
    sci.header.set('EXPTIME', 300.0)

    ad = astrodata.create(phu, [sci])
    ad.filename = f'00000000T000000Z_DDDDE_{arm[0].lower()}_0300.fits'
    return ad


@pytest.fixture(params=['RED', 'BLUE'])
def ad_echelle(request):
    """MaroonX AstroData with synthetic extracted stripes for echelle tests.

    Builds STRIPES, F_STRIPES, and STRIPES_MASKS dicts of sparse matrices
    using real fiber/order positions from the SID lookup files.  Each stripe
    is a flat horizontal band (constant polynomial) so the box-extracted
    flux per column equals ``signal * 2 * slit_height``.
    """
    from pathlib import Path
    from scipy import sparse
    from astropy.table import Table

    arm = request.param
    slit_height = 10
    signal = 100.0
    flat_signal = 1000.0

    # Build the base AstroData object
    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('ORIGNAME', f'00000000T000000Z_DDDDE_{arm[0].lower()}_0300.fits')
    phu.header.set('FIBER1', 'Dark')
    phu.header.set('FIBER2', 'Dark')
    phu.header.set('FIBER3', 'Dark')
    phu.header.set('FIBER4', 'Dark')
    phu.header.set('FIBER5', 'Etalon')

    data = np.full((4400, 4400), signal, dtype=np.float32)
    sci = fits.ImageHDU(data=data, name='SCI')
    sci.header.set('ARM', arm)
    sci.header.set('EXPTIME', 300.0)

    ad = astrodata.create(phu, [sci])
    ad.filename = f'00000000T000000Z_DDDDE_{arm[0].lower()}_0300.fits'
    ad[0].mask = np.zeros((4400, 4400), dtype=np.uint16)

    # Load SID for this arm
    sid_name = 'SID_r.fits' if arm == 'RED' else 'SID_b.fits'
    sid_path = Path(__file__).parent.parent / 'lookups' / 'SID' / sid_name
    sid = Table.read(str(sid_path), hdu='SID')

    flat_data = np.full((4400, 4400), flat_signal, dtype=np.float32)
    mask_data = np.ones((4400, 4400), dtype=int)

    def _extract_stripe(img, poly, sh):
        ny, nx = img.shape
        xx = np.arange(nx)
        y = np.poly1d(poly)(xx)
        sy = np.arange(-sh, sh).repeat(nx).reshape((2 * sh, nx))
        sx = np.tile(np.arange(nx), 2 * sh).reshape((2 * sh, nx))
        idx = np.rint(sy + y).astype(int)
        valid = (idx > 0) & (idx < ny)
        return sparse.coo_matrix(
            (img[idx[valid], sx[valid]], (idx[valid], sx[valid])),
            shape=(ny, nx),
        ).tocsc()

    stripes = {}
    f_stripes = {}
    stripes_masks = {}

    for fiber in sorted(set(sid['identify_fiber'])):
        fkey = f'fiber_{fiber}'
        stripes[fkey] = {}
        f_stripes[fkey] = {}
        stripes_masks[fkey] = {}
        for row in sid[sid['identify_fiber'] == fiber]:
            order = str(row['fiber_order'])
            poly = np.array([float(row['fiber_position'])])
            stripes[fkey][order] = _extract_stripe(data, poly, slit_height)
            f_stripes[fkey][order] = _extract_stripe(flat_data, poly, slit_height)
            stripes_masks[fkey][order] = _extract_stripe(mask_data, poly, slit_height)

    ad[0].STRIPES = stripes
    ad[0].F_STRIPES = f_stripes
    ad[0].STRIPES_MASKS = stripes_masks

    return ad


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
    not_close = ~np.isclose(x, y, rtol=rtol, atol=atol, equal_nan=True)
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
    'FLAT': [
        'N20250701M6126.fits',
        'N20250701M6143.fits',
        'N20250701M6154.fits',
        'N20250701M6164.fits',
        'N20250701M6175.fits',
        'N20250701M6185.fits',
        'N20250701M6215.fits',
        'N20250701M6229.fits',
        'N20250701M6240.fits',
        'N20250701M6250.fits',
        'N20250701M6260.fits',
        'N20250701M6271.fits',
    ],
    'DARK': [
        'N20250717M6066.fits',
        'N20250717M6072.fits',
        'N20250717M6081.fits',
        'N20250717M6091.fits',
        'N20250721M5904.fits',
        'N20250721M5914.fits',
        'N20250721M5930.fits',
        'N20250721M5953.fits',
        'N20250721M5975.fits',
        'N20250721M5997.fits',
        'N20250721M6020.fits',
        'N20250721M6042.fits',
        'N20250721M6059.fits',
        'N20250721M6075.fits',
        'N20250721M6092.fits',
        'N20250721M6108.fits',
        'N20250721M6125.fits',
        'N20250721M6165.fits',
        'N20250721M6206.fits',
        'N20250721M6246.fits',
        'N20250721M6287.fits',
        'N20250721M6328.fits',
        'N20250721M6398.fits',
        'N20250721M6468.fits',
        'N20250721M6539.fits',
        'N20250721M6609.fits',
        'N20250721M6680.fits',
        'N20250721M6780.fits',
        'N20250721M6881.fits',
        'N20250721M6981.fits',
        'N20250721M7081.fits',
        'N20250721M7182.fits',
        'N20250721M7312.fits',
        'N20250721M7443.fits',
        'N20250721M7573.fits',
        'N20250721M7704.fits',
        'N20250721M7835.fits',
        'N20250721M8025.fits',
        'N20250721M8215.fits',
        'N20250721M8406.fits',
        'N20250721M8596.fits',
    ],
    'WAVECAL': [
        'N20250717M5948.fits',
    ],
    'SCIENCE': [
        'N20250717M5299.fits',
    ],
}


def download_mx_file(filename, sub_path='raw_files'):
    """Download a MaroonX raw file from the Gemini Archive.

    Returns the local path on success, or None if the archive responds with
    HTTP 403 (Access Forbidden), typically for proprietary data.
    """
    try:
        return download_from_archive(
            filename, sub_path=sub_path, env_var='DRAGONS_TEST'
        )
    except HTTPError as e:
        warnings.warn(f'{filename}: {e}')
        return None


def download_all_raw_files():
    """Download every file in MAROONX_TEST_MANIFEST into $DRAGONS_TEST/raw_files/."""
    results = {}
    for kind, files in MAROONX_TEST_MANIFEST.items():
        for filename in files:
            results[filename] = download_mx_file(filename)
    return results
