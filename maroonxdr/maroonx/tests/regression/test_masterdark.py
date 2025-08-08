import re
import numpy as np
import pytest
from copy import deepcopy
from pathlib import Path

from astropy.io import fits
import h5py

import astrodata
from gempy.adlibrary import dataselect

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


# =========================================================
# TESTS
# =========================================================

def test_masterdark(legacy_darks_path, science_dir):

    old_file = legacy_darks_path / "20241115T19_masterdark_mean_DDDDE_b_0300.fits"
    
    # get all flat files
    raw_files = sorted([str(f) for f in science_dir.glob('*.fits')])
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


@pytest.mark.slow
def test_fitDarkCoefficients(arm, legacy_darks_path, processed_dark_path):

    legacy_file = legacy_darks_path / f"masterdarks_coeffs_202411xx_{arm.lower()}.npz"
    legacy_coeffs = np.load(legacy_file)

    # get all dark files
    all_files = sorted([str(f) for f in processed_dark_path.glob('*.fits')])
    selected_dark = dataselect.select_data(all_files, tags=['PROCESSED', 'DARK', arm])

    # read files and instantiate the primitive class
    adinput = [astrodata.open(f) for f in selected_dark]
    
    # copy the file for trimming
    # need to copy a raw file that does not have the trim overscan keyword applied
    # legacy_z0_ad = deepcopy(adinput[0])
    # legacy_z1_ad = deepcopy(adinput[0])
    # legacy_z0_ad[0].data = legacy_coeffs['z0']
    # legacy_z1_ad[0].data = legacy_coeffs['z1']

    p = MAROONX(adinput)
    p.prepare()
    p.checkArm()
    p.checkMaster()
    # p.checkND()
    adout = p.fitDarkCoefficients()

    z0 = adout[0][0].COEFF_Z0
    z1 = adout[0][0].COEFF_Z1
    logexptime = adout[0][0].LOGEXPTIME

    # Trim legacy data
    legacy_z0 = MAROONX([legacy_z0_ad]).trimOverscan()[0][0].data
    legacy_z1 = MAROONX([legacy_z1_ad]).trimOverscan()[0][0].data

    # Compare the data
    # np.testing.assert_allclose(z0, legacy_z0, rtol=1e-5, atol=1e-5, 
    #     err_msg='Data mismatch')
    # np.testing.assert_allclose(z1, legacy_z1, rtol=1e-5, atol=1e-5, 
    #     err_msg='Data mismatch')
    
    np.testing.assert_allclose(logexptime['logexptime'].value, legacy_coeffs['logexptime'])



# =====================================================
# Helper functions
# =====================================================


def old_trim():
    # old dark coefficient need trimminig before comparisson
    ...