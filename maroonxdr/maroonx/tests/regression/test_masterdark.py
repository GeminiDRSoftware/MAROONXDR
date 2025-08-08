import re
import numpy as np
import pytest
from pathlib import Path

from astropy.io import fits
import h5py

import astrodata
from gempy.adlibrary import dataselect

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


# =========================================================
# FILES TO COMPARE
# =========================================================

# Paths to masterdarks coefficient files
# COEFF_OUTPUT_PATH = Path("/home/martin/Documentos/Projects/MaroonX/maroonx_base/data2/MaroonX_spectra_reduced/Maroonx_masterframes/202411xx/darks/")
# COEFF_REFERENCE_PATH = Path("/home/martin/Documentos/Projects/MaroonX/resources/Martín_setup/darks/")

# Paths to real masterdarks files
OLD_FILES_PATH = Path("/home/martin/Projects/MaroonX/legacy/maroonx_base/data2/MaroonX_spectra_reduced/Maroonx_masterframes/202411xx/darks")
NEW_FILES_PATH = Path("/home/martin/Projects/MaroonX/MAROONXDR/calibrations/processed_dark")

SCIENCE_DIR = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir')

PROCESSED_DARK = Path('/home/martin/Projects/MaroonX/MAROONXDR/calibrations/processed_dark')

# =========================================================
# TESTS
# =========================================================

def test_masterdark():

    old_file = OLD_FILES_PATH / "20241115T19_masterdark_mean_DDDDE_b_0300.fits"
    
    # get all flat files
    raw_files = sorted([str(f) for f in SCIENCE_DIR.glob('*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'DARK', 'BLUE', '300s'])

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


    with fits.open(old_file) as old_hdu:
        old_data = old_hdu[0].data
        new_data = adouts[0][0].data

        assert old_data.shape == new_data.shape, f"Shape mismatch: {old_data.shape} != {new_data.shape}"

        # Compare the data
        np.testing.assert_allclose(old_data, new_data, rtol=1e-5, atol=1e-5, 
            err_msg='Data mismatch')


@pytest.mark.parametrize("arm", ["BLUE", "RED"])
def test_fitDarkCoefficients(arm):

    old_file = OLD_FILES_PATH / f"masterdarks_coeffs_202411xx_{arm.lower()}.npz"
    old_coeffs = np.load(old_file)

    # get all dark files
    all_files = sorted([str(f) for f in PROCESSED_DARK.glob('*.fits')])
    selected_dark = dataselect.select_data(all_files, tags=['PROCESSED', 'DARK', 'BLUE'])

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

    # Compare the data
    # np.testing.assert_allclose(z0, old_coeffs['z0'], rtol=1e-5, atol=1e-5, 
    #     err_msg='Data mismatch')
    # np.testing.assert_allclose(z1, old_coeffs['z1'], rtol=1e-5, atol=1e-5, 
    #     err_msg='Data mismatch')
    
    old_coeffs['logexptime']
    logexptime['logexptime'].value
    np.testing.assert_allclose(logexptime['logexptime'].value, old_coeffs['logexptime'])



# =====================================================
# Helper functions
# =====================================================


def old_trim():
    # old dark coefficient need trimminig before comparisson
    ...