

from copy import deepcopy
import numpy as np
import pandas as pd
import pytest

from scipy import sparse
from pathlib import Path

import astrodata
from astropy.io import fits
import h5py

from gempy.adlibrary import dataselect
from gempy.utils import logutils

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


# =========================================================
# FILES TO COMPARE
# =========================================================

OLD_BASE_PATH = Path("/home/martin/Projects/MaroonX/legacy/maroonx_base/data2/")
OLD_FILES_PATH = OLD_BASE_PATH / Path("MaroonX_spectra_reduced/20241124")
OLD_FLAT_FILES_PATH = OLD_BASE_PATH / Path("MaroonX_spectra_reduced/Maroonx_masterframes/202411xx/flats")

SCIENCE_DIR = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir')

# Set logger
logutils.config(file_name="test_extractions.log", mode="debug", stomp=True)
log = logutils.get_logger("test_extractions")
log.setLevel("DEBUG")

# =========================================================
# TESTS
# =========================================================
USE_CACHE = True

SCIENCE_FILES = [
    '20241124T041907Z_SOOOE_b_0300',
]


@pytest.mark.parametrize("arm", ["BLUE"])
def test_extractStripes(arm):

    old_file = OLD_FILES_PATH / "20241124T162336Z_DEEEE_b_0030.hdf"

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in SCIENCE_DIR.glob('20241124T162336Z_DEEEE_*.fits')])
    
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'WAVECAL', arm])

    # Primitives
    adinput = [astrodata.open(f) for f in selected_spect]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)
    
    adout = p.extractStripes()

    # STRIPES extension is an intermediate dictionary of fibers and orders
    new_stripes = adout[0][0].STRIPES

    fail_counter = 0
    for f in new_stripes.keys():
        for o in new_stripes[f].keys():
            # this function comes from the legacy pipeline
            mat = load_sparse_mat(f'extracted_stripes/{f}/{o}', str(old_file))
            
            legacy_stripe = mat.toarray()

            new_stripe = new_stripes[f][o].toarray()
            
            try:
                np.testing.assert_allclose(legacy_stripe, new_stripe)
                print(f'fiber/order : {f}/{o} [OK]')
            except AssertionError as err:
                fail_counter += 1
                print(f'fiber/order : {f}/{o} [FAIL]')
    assert fail_counter == 0


@pytest.mark.parametrize("arm", ["BLUE"])
def test_extractStripes_fromFlat(arm):

    old_flat = OLD_FLAT_FILES_PATH / "20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf"

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in SCIENCE_DIR.glob('20241124T162336Z_DEEEE_*.fits')])
    
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'WAVECAL', arm])

    # Primitives
    adinput = [astrodata.open(f) for f in selected_spect]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)
    
    adout = p.extractStripes()

    # STRIPES extension is an intermediate dictionary of fibers and orders    
    
    new_flat_stripes = adout[0][0].F_STRIPES

    fail_counter = 0
    for f in new_flat_stripes.keys():
        for o in new_flat_stripes[f].keys():
            # this function comes from the legacy pipeline
            flat_mat = load_sparse_mat(f'extracted_stripes/{f}/{o}', str(old_flat))
            
            legacy_flat_stripe = flat_mat.toarray()

            new_flat_stripe = new_flat_stripes[f][o].toarray()
            
            try:
                np.testing.assert_allclose(legacy_flat_stripe, new_flat_stripe)
                print(f'fiber/order : {f}/{o} [OK]')
            except AssertionError as err:
                fail_counter += 1
                print(f'fiber/order : {f}/{o} [FAIL]')
    assert fail_counter == 0


@pytest.mark.parametrize("arm", ["BLUE"])
def test_extractStripes_fromScience(arm):

    # old_flat = OLD_FLAT_FILES_PATH / "20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf"
    
    old_science = OLD_FILES_PATH / "20241124T041907Z_SOOOE_b_0300.hdf"
    LEGACY_TEST_NPY_PATH = OLD_BASE_PATH.parent / "legacy_bkg_arrays"
    legacy_npy = LEGACY_TEST_NPY_PATH / "20241124T041907Z_SOOOE_b_0300_test.npy"
    legacy_dict = np.load(legacy_npy, allow_pickle=True).item()
    # legacy_dict.keys() has the following keys:
    # dict_keys(['raw', 'bias_corrected', 'overscan_removed', 'orientation_corrected', 'dark_science', 'back_var', 'index_fiber', 'index_order', 'mask', 'science_straylight_removed'])

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in SCIENCE_DIR.glob('20241124T041907Z_SOOOE_*.fits')])
    
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'SCI', arm])

    # Primitives
    adinput = [astrodata.open(f) for f in selected_spect]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()

    # p.addVAR() --> SYNTH DARK VARIANCE ?

    # test raw data is the same as legacy
    np.testing.assert_allclose(
        p.streams["main"][0][0].data, 
        legacy_dict['raw'], 
        rtol=1e-4, atol=1e-4)

    p.subtractOverscan()
    # test overscan removed data is the same as legacy
    np.testing.assert_allclose(
        p.streams["main"][0][0].data, 
        legacy_dict['bias_corrected'], 
        rtol=1e-4, atol=1e-4)

    p.trimOverscan()
    # test overscan removed data is the same as legacy
    np.testing.assert_allclose(
        p.streams["main"][0][0].data, 
        legacy_dict['overscan_removed'], 
        rtol=1e-4, atol=1e-4)

    p.correctImageOrientation()
    # test orientation corrected data is the same as legacy
    np.testing.assert_allclose(
        p.streams["main"][0][0].data, 
        legacy_dict['orientation_corrected'], 
        rtol=1e-4, atol=1e-4)

    p.removeStraylight()
    # test straylight corrected data is the same as legacy
    np.testing.assert_allclose(
        p.streams["main"][0][0].data, 
        legacy_dict['science_straylight_removed'], 
        rtol=1e-4, atol=1e-4)
    # ====================================================== OK line

    p.extractStripes()

    # STRIPES extension is an intermediate dictionary of fibers and orders
    new_stripes = p.streams["main"][0][0].STRIPES

    fail_counter = 0
    for f in new_stripes.keys():
        for o in new_stripes[f].keys():
            # this function comes from the legacy pipeline
            mat = load_sparse_mat(f'extracted_stripes/{f}/{o}', str(old_science))
            
            legacy_stripe = mat.toarray()

            new_stripe = new_stripes[f][o].toarray()
            
            try:
                np.testing.assert_allclose(legacy_stripe, new_stripe)
                print(f'fiber/order : {f}/{o} [OK]')
            except AssertionError as err:
                fail_counter += 1
                print(f'fiber/order : {f}/{o} [FAIL]')
    assert fail_counter == 0

    # ==============================================================



@pytest.mark.parametrize("arm", ["BLUE"])
def test_boxExtraction(arm):

    old_file = OLD_FILES_PATH / "20241124T162336Z_DEEEE_b_0030.hdf"

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in SCIENCE_DIR.glob('20241124T162336Z_DEEEE_*.fits')])
    
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'WAVECAL', arm])

    # Primitives
    adinput = [astrodata.open(f) for f in selected_spect]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)
    
    p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
    
    adout = p.boxExtraction() # extracts spectra from stripes

    legacy_box = load_dict_from_hdf5(str(old_file), "box_extraction/")

    fail_counter = 0
    for fiber in range(1, 6):
        new_box = getattr(adout[0][0], f"BOX_REDUCED_FIBER_{fiber}")
        if new_box.size == 1:
            # only fiber 1 should be empty
            assert fiber == 1, "Only fiber 1 should be empty in box extraction"
            continue

        orders = getattr(adout[0][0], f"REDUCED_ORDERS_FIBER_{fiber}")
        orders = orders.astype(int)

        for idx, order in enumerate(orders):
            legacy_order = legacy_box[f"fiber_{fiber}"][f"{order}"]
            new_order = new_box[idx, :]
            #np.testing.assert_allclose(legacy_order, new_order, rtol=1e-4)

            try:
                np.testing.assert_allclose(legacy_order, new_order, rtol=1e-4)
                # print(f'fiber/order : {fiber}/{order} [OK]')
            except AssertionError as err:
                fail_counter += 1
                # print(f'fiber/order : {fiber}/{order} [FAIL]')
    assert fail_counter == 0

@pytest.mark.parametrize("science_filename", SCIENCE_FILES)
def test_optimalExtraction(science_filename):
    """
    This test is not implemented yet, as the optimal extraction is not yet implemented in MaroonX.
    """

    old_file = OLD_FILES_PATH / (science_filename + ".hdf")

    if USE_CACHE:
        # Use previously saved data on science_dir
        adout = [astrodata.open(SCIENCE_DIR / (science_filename + "_reduced.fits"))]
    else:
        # Primitives
        adinput = [astrodata.open(SCIENCE_DIR / (science_filename + ".fits"))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        # get and save wavelength solution (either static reference or frame's unique sim cal solved)
        p.darkSubtraction()
        p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra TODO Skip dark for fiber 5
        
        adout = p.optimalExtraction() # extracts spectra from stripes

    legacy_box = load_dict_from_hdf5(str(old_file), "box_extraction/")
    legacy_opt = load_dict_from_hdf5(str(old_file), "optimal_extraction/")
    legacy_var = load_dict_from_hdf5(str(old_file), "optimal_var/")

    fail_counter = 0
    for fiber in range(1, 6):
        new_box = getattr(adout[0][0], f"BOX_REDUCED_FIBER_{fiber}")
        new_opt = getattr(adout[0][0], f"OPTIMAL_REDUCED_FIBER_{fiber}")
        if new_opt.size == 1:
            assert fiber in [1, 5], "Only fibers 1 and 5 should be empty in optimal extraction"
            continue

        orders = getattr(adout[0][0], f"REDUCED_ORDERS_FIBER_{fiber}")
        orders = orders.astype(int)

        for idx, order in enumerate(orders):
            legacy_order = legacy_opt[f"fiber_{fiber}"][f"{order}"]
            new_order = new_opt[idx, :]
            #np.testing.assert_allclose(legacy_order, new_order, rtol=1e-4)
            try:
                np.testing.assert_allclose(legacy_order, new_order, rtol=1e-4)
                print(f'fiber/order : {fiber}/{order} [OK]')
            except AssertionError as err:
                fail_counter += 1
                print(f'fiber/order : {fiber}/{order} [FAIL]')

            legacy_order = legacy_box[f"fiber_{fiber}"][f"{order}"]
            new_order = new_box[idx, :]
            try:
                np.testing.assert_allclose(legacy_order, new_order, rtol=1e-4)
                print(f'fiber/order [box]: {fiber}/{order} [OK]')
            except AssertionError as err:
                fail_counter += 1
                print(f'fiber/order [box]: {fiber}/{order} [FAIL]')
    assert fail_counter == 0


@pytest.mark.parametrize("arm", ["BLUE"])
def test_staticWavelengthSolution(arm):

    old_file = OLD_FILES_PATH / "20241124T162336Z_DEEEE_b_0030.hdf"

    # Load old peak data. columns are lowercase
    legacy_wls = load_dict_from_hdf5(str(old_file), 'wavelengths_static/')

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in SCIENCE_DIR.glob('20241124T162336Z_DEEEE_*.fits')])
    
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'WAVECAL', arm])

    # Primitives
    adinput = [astrodata.open(f) for f in selected_spect]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)
    
    p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
    p.boxExtraction() # extracts spectra from stripes

    #p.getPeaksAndPolynomials() # not eeded
    adout = p.staticWavelengthSolution()

    fail_counter = 0
    for fiber in range(1, 6):
        new_wls = getattr(adout[0][0], f"WLS_STATIC_FIBER_{fiber}")
        if new_wls.size == 1:
            # non reduced fibers (usually fiber 1) have size 1
            continue

        orders = getattr(adout[0][0], f"REDUCED_ORDERS_FIBER_{fiber}")
        orders = orders.astype(int)

        for idx, order in enumerate(orders):
            legacy_order_wls = legacy_wls[f"fiber_{fiber}"][f"{order}"]
            new_order_wls = new_wls[idx, :]

            try:
                np.testing.assert_allclose(legacy_order_wls, new_order_wls, rtol=1e-4)
                # print(f'fiber/order : {fiber}/{order} [OK]')
            except AssertionError as err:
                fail_counter += 1
                # print(f'fiber/order : {fiber}/{order} [FAIL]')
    assert fail_counter == 0

 


# =====================================================
# HDF5 helper functions
# =====================================================

def load_mat(name, filepath):
    with h5py.File(filepath, mode="r") as h5f:
        return h5f[name][()]

def load_dict_from_hdf5(filename, path):
    if not path.endswith("/"):
        path += "/"
    if isinstance(filename, str):
        with h5py.File(filename, 'r', libver='latest') as h5file:
            return recursively_load_dict_contents_from_group(h5file, path)
    else:
        return recursively_load_dict_contents_from_group(filename, path)

def recursively_load_dict_contents_from_group(h5file, path):
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans

def load_sparse_mat(name, store='store.h5'):
    with h5py.File(store, mode="r") as f:
        pars = []
        for par in ('data', 'indices', 'indptr'):
            pars.append(f['%s/%s' % (name, par)][()])
        pars.append(f[name].attrs['h5sparse_shape'])
    m = sparse.csc_matrix(tuple(pars[:3]), shape=pars[3])
    return m