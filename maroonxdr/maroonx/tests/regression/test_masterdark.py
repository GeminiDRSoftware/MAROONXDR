import re
import numpy as np
import pytest
from copy import deepcopy
from pathlib import Path

from astropy.io import fits
import h5py

import astrodata
from gempy.adlibrary import dataselect
from gempy.utils import logutils

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX
#from maroonxdr.maroonx.primitives_maroonx_echelle import create_synthetic_dark

# Set logger
logutils.config(file_name="test_masterdarks.log", mode="debug", stomp=True)
log = logutils.get_logger("test_masterdarks.log")
log.setLevel("DEBUG")

# =========================================================
# TESTS
# =========================================================

@pytest.mark.parametrize("exptime", ["60", "120", "300", "600", "900", "1200", "1800"])
def test_master_dark(arm, legacy_darks_path, exptime):

    legacy_suffix =  f"masterdark_mean_DDDDE_{arm[0].lower()}_{exptime.zfill(4)}"
    legacy_file = list(legacy_darks_path.glob(f"202411*T*_{legacy_suffix}.fits"))
    assert len(legacy_file) == 1, f"Expected 1 file, got {len(legacy_file)}"
    legacy_file = list(legacy_file)[0]

    # get all flat files
    raw_files = sorted([str(f) for f in Path().glob('*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'DARK', arm, f'{exptime}s'])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.subtractOverscan()
    # No trim overscan and no correct image orientation
    p.addVAR(read_noise=True, poisson_noise=True)
    adouts = p.stackDarks(scale_mode='first_frame', lsigma=2.0, hsigma=2.0)


    with fits.open(legacy_file) as old_hdu:
        old_data = old_hdu[0].data
        new_data = adouts[0][0].data

        assert old_data.shape == new_data.shape, f"Shape mismatch: {old_data.shape} != {new_data.shape}"

        # Compare the data
        np.testing.assert_allclose(old_data, new_data, rtol=1e-5, atol=1e-5, 
            err_msg='Data mismatch')


@pytest.mark.parametrize("exptime", ["60", "120", "300", "600", "900", "1200", "1800"])
def test_master_dark_with_legacy_post_trim(arm, legacy_darks_path, exptime, ad_empty_dark):

    legacy_suffix =  f"masterdark_mean_DDDDE_{arm[0].lower()}_{exptime.zfill(4)}"
    legacy_file = list(legacy_darks_path.glob(f"202411*T*_{legacy_suffix}.fits"))
    assert len(legacy_file) == 1, f"Expected 1 file, got {len(legacy_file)}"
    legacy_file = list(legacy_file)[0]

    # get all flat files
    raw_files = sorted([str(f) for f in Path().glob('*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'DARK', arm, f'{exptime}s'])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_spect]
    p = MAROONX(adinput)

    p.prepare()
    p.checkArm()
    p.checkND()
    p.addDQ()
    p.subtractOverscan()
    # Trimming and flipping
    p.trimOverscan()
    p.correctImageOrientation()

    p.addVAR(read_noise=True, poisson_noise=True)
    adouts = p.stackDarks(scale_mode='first_frame', lsigma=2.0, hsigma=2.0)

    with fits.open(legacy_file) as old_hdu:
        old_data = old_hdu[0].data
        new_data = adouts[0][0].data

        ad_empty_dark[0].data = old_data.copy()     # not empty anymore
        ad_dark = trim_and_flip_legacy(ad_empty_dark)

        old_data_post_trim = ad_dark[0].data
        assert old_data_post_trim.shape == new_data.shape, f"Shape mismatch: {old_data_post_trim.shape} != {new_data.shape}"

        # Compare the data
        np.testing.assert_allclose(old_data_post_trim, new_data, rtol=1e-2, atol=1e-2, 
            err_msg='Data mismatch')

@pytest.mark.slow
def test_fitDarkCoefficients(arm, legacy_darks_path, processed_dark_path, ad_empty_dark):

    legacy_file = legacy_darks_path / f"masterdarks_coeffs_202411xx_{arm.lower()}.npz"
    legacy_coeffs = np.load(legacy_file)

    # get all dark files
    all_files = sorted([str(f) for f in processed_dark_path.glob('*.fits')])
    selected_dark = dataselect.select_data(all_files, tags=['PROCESSED', 'DARK', arm])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_dark]

    p = MAROONX(adinput)
    p.prepare()
    p.checkArm()
    p.checkMaster()
    # p.checkND()
    adout = p.fitDarkCoefficients()

    z0 = adout[0][0].COEFF_Z0
    z1 = adout[0][0].COEFF_Z1
    logexptime = adout[0][0].LOGEXPTIME

    # Load legacy data into an empty astrodata object and trim it
    # legacy_z0_ad = deepcopy(ad_empty_dark)
    # legacy_z1_ad = deepcopy(ad_empty_dark)
    # legacy_z0_ad[0].data = legacy_coeffs['z0']
    # legacy_z1_ad[0].data = legacy_coeffs['z1']
    
    # legacy_z0_p = MAROONX([legacy_z0_ad])
    # legacy_z1_p = MAROONX([legacy_z1_ad])
    
    # legacy_z0_p.trimOverscan()
    # legacy_z1_p.trimOverscan()
    
    # legacy_z0_p.correctImageOrientation()
    # legacy_z1_p.correctImageOrientation()

    # legacy_z0 = legacy_z0_p.streams['main'][0][0].data
    # legacy_z1 = legacy_z1_p.streams['main'][0][0].data
    
    legacy_z0 = legacy_coeffs['z0']
    legacy_z1 = legacy_coeffs['z1']

    # Compare the data
    np.testing.assert_allclose(z0, legacy_z0, rtol=1e-2, atol=1e-2, err_msg='Data mismatch')
    np.testing.assert_allclose(z1, legacy_z1, rtol=1e-2, atol=1e-2, err_msg='Data mismatch')
    np.testing.assert_allclose(logexptime['logexptime'].value, legacy_coeffs['logexptime'])


#@pytest.mark.parametrize("exptime", ["300"])
def test_create_synthetic_masterdark(arm, legacy_darks_path, processed_dark_path):
    """
    Test the creation of a synthetic master dark.
    """

    exptime = "300"
    ndfilter = 119.05160812

    legacy_suffix =  f"masterdark_mean_DDDDE_{arm[0].lower()}_{exptime.zfill(4)}"
    legacy_file = list(legacy_darks_path.glob(f"202411xx_{legacy_suffix}.fits"))
    assert len(legacy_file) == 1, f"Expected 1 file, got {len(legacy_file)}"
    legacy_file = legacy_file[0]


    # get coeff dark files
    all_files = sorted([str(f) for f in processed_dark_path.glob('*.fits')])
    selected_dark = dataselect.select_data(all_files, tags=['DARK_COEFF', arm])
    assert len(selected_dark) == 1, "Expected 1 file with tag DARK_COEFF"

    ad_coeff = astrodata.open(selected_dark[0])

    # create synthetic master dark
    synthetic_dark, actual_exptime, actual_nd, factor = create_synthetic_dark(ad_coeff, float(exptime), ndfilter)
    print(actual_exptime, actual_nd, factor)

    with fits.open(legacy_file) as hdul:
        legacy_dark = hdul[0].data

        # check that the synthetic dark frame is within 5% of the legacy dark frame
        np.testing.assert_allclose(synthetic_dark, legacy_dark, rtol=1e-5)

# =========================================================================================
def trim_and_flip_legacy(ad):
    """
    Trim the input AstroData object to a 2D array and flip it vertically.
    """
    # Trim the input AstroData object to a 2D array
    p = MAROONX([ad])
    p.trimOverscan()
    p.correctImageOrientation()
    return p.streams['main'][0].data

def create_synthetic_dark(ad_coeff, exptime_value=None, nd_value=None):
    """
    Create a synthetic dark frame from coefficients file.
    
    This function generates a synthetic dark frame using pre-computed coefficients
    and the log-linear relationship: dark = z1 + z0 * log10(exptime * factor)
    
    Parameters
    ----------
    ad_coeff : AstroData
        AstroData object containing coefficient data with extensions:
        - COEFF_Z0: slope coefficients for each pixel
        - COEFF_Z1: intercept coefficients for each pixel  
        - LOGEXPTIME: table with exposure times and ND filter data
    exptime_value : float, optional
        Target exposure time in seconds. If None, must provide nd_value.
    nd_value : float, optional
        Target ND filter position. If None, will be calculated from exptime_value.
        
    Returns
    -------
    synthetic_dark : numpy.ndarray
        2D array containing the synthetic dark frame
    actual_exptime : float
        The actual exposure time used (may differ from input due to ND corrections)
    actual_nd : float
        The actual ND filter position used
    factor : float
        The correction factor applied (1.0 if no correction needed)
        
    Raises
    ------
    ValueError
        If neither exptime_value nor nd_value is provided, or if required 
        extensions are missing from ad_coeff
    """
    
    # Validate inputs
    if exptime_value is None and nd_value is None:
        raise ValueError("Either exptime_value or nd_value must be provided")
    
    # Check for required extensions
    required_extensions = ['COEFF_Z0', 'COEFF_Z1', 'LOGEXPTIME']
    for ext_name in required_extensions:
        if not hasattr(ad_coeff[0], ext_name):
            raise ValueError(f"Required extension {ext_name} not found")
    
    # Extract coefficient arrays
    z0 = ad_coeff[0].COEFF_Z0
    z1 = ad_coeff[0].COEFF_Z1
    logexptime_table = ad_coeff[0].LOGEXPTIME
    
    # Extract calibration data
    logexptimes = np.array(logexptime_table['logexptime'])
    exptimes = np.array(logexptime_table['exptime'])
    
    # Check if ND filter data is available
    has_nd_data = 'ndfilter' in logexptime_table.colnames
    if has_nd_data:
        ndfilters = np.array(logexptime_table['ndfilter'])
    else:
        ndfilters = np.zeros_like(exptimes)  # Default to zero if no ND data
    
    # Initialize variables
    factor = 1.0
    actual_exptime = exptime_value
    actual_nd = nd_value
    
    # Case 1: Only exposure time provided
    if exptime_value is not None and nd_value is None:
        actual_exptime = exptime_value
        
        if has_nd_data and len(np.unique(ndfilters)) > 1:
            # Calculate expected ND position from exposure time
            z_nd_logt = np.polyfit(logexptimes, ndfilters, 1)
            f_nd_logt = np.poly1d(z_nd_logt)
            actual_nd = f_nd_logt(np.log10(exptime_value))
        else:
            actual_nd = ndfilters[0] if len(ndfilters) > 0 else 0.0
    
    # Case 2: Only ND value provided  
    elif nd_value is not None and exptime_value is None:
        actual_nd = nd_value
        
        if has_nd_data and len(np.unique(ndfilters)) > 1:
            # Calculate exposure time from ND position
            z_logt_nd = np.polyfit(ndfilters, logexptimes, 1)
            f_logt_nd = np.poly1d(z_logt_nd)
            actual_exptime = 10**(f_logt_nd(nd_value))
        else:
            # Use median exposure time as default
            actual_exptime = np.median(exptimes)
    
    # Case 3: Both exposure time and ND value provided
    elif exptime_value is not None and nd_value is not None:
        actual_exptime = exptime_value
        actual_nd = nd_value
        
        # Check for ND filter mismatch and apply correction if needed
        if has_nd_data and len(np.unique(ndfilters)) > 1:
            z_nd_logt = np.polyfit(logexptimes, ndfilters, 1)
            f_nd_logt = np.poly1d(z_nd_logt)
            z_logt_nd = np.polyfit(ndfilters, logexptimes, 1) 
            f_logt_nd = np.poly1d(z_logt_nd)
            
            expected_nd = f_nd_logt(np.log10(exptime_value))
            nd_difference = abs(actual_nd - expected_nd)
            
            # Apply correction factor if ND mismatch is significant (>0.2)
            if nd_difference > 0.2:
                logt_nominal = f_logt_nd(actual_nd)
                factor = 10**(np.log10(exptime_value) - logt_nominal)
    
    # Calculate synthetic dark frame
    # Formula: dark = z1 + z0 * log10(exptime * factor)
    effective_logexptime = np.log10(actual_exptime * factor)
    synthetic_dark = z1 + z0 * effective_logexptime
    
    return synthetic_dark.astype(np.float32), actual_exptime, actual_nd, factor