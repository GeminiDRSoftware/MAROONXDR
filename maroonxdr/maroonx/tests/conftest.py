import os
import pytest
from pathlib import Path
import tempfile
import shutil
from contextlib import contextmanager

import astrodata
from gempy.adlibrary import dataselect
import maroonx_instruments  # noqa : important to load adclass tags

# =========================================================
# DRAGONS TEST CONFIGURATION
# =========================================================

def get_maroonx_dragons_test_path():
    """
    Get the DRAGONS_TEST environment variable path.
    
    Returns
    -------
    Path or None
        Path to DRAGONS test data directory, or None if not set
    """
    p = os.environ.get("MAROONX_DRAGONS_TEST")
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
    p = os.environ.get("MAROONX_LEGACY_TEST")
    return Path(p) if p else None


# =========================================================
# PATH FIXTURES
# =========================================================

@pytest.fixture(scope="session")
def dragons_test_root():
    """
    Fixture providing the root DRAGONS test data directory.
    
    Returns
    -------
    Path
        Root path to DRAGONS test data
    """
    root = get_maroonx_dragons_test_path()
    if root is None:
        pytest.skip("MAROONX_DRAGONS_TEST environment variable not set")
    if not root.exists():
        pytest.skip(f"DRAGONS test root does not exist: {root}")
    return root

@pytest.fixture(scope="session")
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
        pytest.skip("MAROONX_LEGACY_TEST environment variable not set")
    if not root.exists():
        pytest.skip(f"Legacy root does not exist: {root}")
    return root

@pytest.fixture(scope="function")
def legacy_darks_path(legacy_test_root):
    """
    Fixture providing path to test legacy data directory.
    """
    path = legacy_test_root / "MaroonX_spectra_reduced" / "Maroonx_masterframes" / "202411xx" / "darks"
    if not path.exists():
        pytest.skip(f"Legacy data directory does not exist: {path}")
    return path

@pytest.fixture(scope="function")
def legacy_flats_path(legacy_test_root):
    """
    Fixture providing path to test legacy data directory.
    """
    path = legacy_test_root / "MaroonX_spectra_reduced" / "Maroonx_masterframes" / "202411xx" / "flats"
    if not path.exists():
        pytest.skip(f"Legacy data directory does not exist: {path}")
    return path

@pytest.fixture(scope="function")
def science_dir(dragons_test_root):
    """
    Fixture providing path to test legacy data directory.
    """
    path = dragons_test_root / "science_dir"
    if not path.exists():
        pytest.skip(f"Science dir directory does not exist: {path}")
    return path

@pytest.fixture(scope="function")
def processed_dark_path(dragons_test_root):
    """
    Fixture providing path to test processed dark calibrations.
    """
    path = dragons_test_root / "calibrations" / "processed_dark"
    if not path.exists():
        pytest.skip(f"Science dir directory does not exist: {path}")
    return path

# =========================================================
# FIXTURES
# =========================================================

@pytest.fixture(params=["BLUE", "RED"])
def arm(request):
    """
    Fixture providing color tag name for the arm.
    """    
    return request.param

# =========================================================
# CONFIGURATION HOOKS
# =========================================================

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
        "markers", 
        "regression: marks tests as regression tests (deselect with '-m \"not regression\"')"
    )
    config.addinivalue_line(
        "markers",
        "slow: mark test as slow tests (deselect with '-m \"not slow\"')"
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
        if "regression" in str(item.path):
            item.add_marker(pytest.mark.regression)