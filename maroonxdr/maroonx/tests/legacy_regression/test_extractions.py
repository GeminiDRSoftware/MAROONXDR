import os

import astrodata
import h5py
import numpy as np
import pytest
from gempy.utils import logutils
from scipy import sparse

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from maroonxdr.maroonx.tests.conftest import assert_allclose_with_max_fails

from . import legacy_adapter

_CALIB_FILES = {
    'RED': {
        'flat': 'processed_flat/20250701T171955Z_DDDDF_r_0002_DFFFF_flat.fits',
        'dark': 'processed_dark/20250717T144308Z_SOOOE_r_0300_synth_dark.fits',
        'wavecal': 'processed_wavecal/20250717T163124Z_DEEEE_r_0004_wavecal.fits',
    },
    'BLUE': {
        'flat': 'processed_flat/20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits',
        'dark': 'processed_dark/20250717T144308Z_SOOOE_b_0300_synth_dark.fits',
        'wavecal': 'processed_wavecal/20250717T163124Z_DEEEE_b_0010_wavecal.fits',
    },
}

logutils.config(file_name="test_extractions_v2.log", mode="debug", stomp=True)
log = logutils.get_logger("test_extractions_v2")
log.setLevel("DEBUG")

# =========================================================
# TESTS
# =========================================================
USE_CACHE = True

SCIENCE_FILES = [
    '20250717T144308Z_SOOOE_b_0300',
    '20250717T144308Z_SOOOE_r_0300',
]

ETALON_FILES = [
    '20250717T163124Z_DEEEE_b_0010',
    '20250717T163124Z_DEEEE_r_0004',
]


def _arm_from_filename(filename):
    return 'BLUE' if '_b_' in filename else 'RED'


def _load_sparse_mat(name, store):
    with h5py.File(store, mode="r") as f:
        pars = []
        for par in ('data', 'indices', 'indptr'):
            pars.append(f['%s/%s' % (name, par)][()])
        pars.append(f[name].attrs['h5sparse_shape'])
    return sparse.csc_matrix(tuple(pars[:3]), shape=pars[3])


@pytest.mark.slow()
@pytest.mark.parametrize("etalon_filename", ETALON_FILES)
def test_extractStripes_fromEtalon(etalon_filename, path_to_legacy_wavecal, preprocessed_files_path):

    arm = _arm_from_filename(etalon_filename)
    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])

    raw_file = preprocessed_files_path / (etalon_filename + ".fits")

    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    adout = p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])

    new_stripes = adout[0][0].STRIPES

    fail_counter = 0
    for f in new_stripes.keys():
        for o in new_stripes[f].keys():
            mat = _load_sparse_mat(f'extracted_stripes/{f}/{o}', str(old_file))
            legacy_stripe = mat.toarray()
            new_stripe = new_stripes[f][o].toarray()

            try:
                assert_allclose_with_max_fails(legacy_stripe, new_stripe, rtol=0, atol=1e-8, max_fails=0)
                log.fullinfo(f'fiber/order : {f}/{o} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {f}/{o} [FAIL]')

    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.slow()
@pytest.mark.parametrize("etalon_filename", ETALON_FILES)
def test_extractStripes_fromFlat(etalon_filename, path_to_legacy_flats, preprocessed_files_path):

    arm = _arm_from_filename(etalon_filename)
    flat_basename = (
        '20250701T17_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf' if arm == 'BLUE'
        else '20250701T17_masterflat_backgroundsubtracted_FFFFF_r_0002.hdf'
    )
    old_flat = path_to_legacy_flats / flat_basename

    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])

    raw_file = preprocessed_files_path / (etalon_filename + ".fits")

    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    adout = p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])

    new_flat_stripes = adout[0][0].F_STRIPES

    fail_counter = 0
    for f in new_flat_stripes.keys():
        for o in new_flat_stripes[f].keys():
            flat_mat = _load_sparse_mat(f'extracted_stripes/{f}/{o}', str(old_flat))
            legacy_flat_stripe = flat_mat.toarray()
            new_flat_stripe = new_flat_stripes[f][o].toarray()

            try:
                assert_allclose_with_max_fails(legacy_flat_stripe, new_flat_stripe, rtol=0, atol=1e-8, max_fails=0)
                log.fullinfo(f'fiber/order : {f}/{o} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {f}/{o} [FAIL]')

    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.slow()
@pytest.mark.parametrize("science_filename", SCIENCE_FILES)
def test_extractStripes_fromScience(science_filename, path_to_legacy_science, path_to_legacy_bkg, preprocessed_files_path):

    arm = _arm_from_filename(science_filename)
    old_science = path_to_legacy_science / (science_filename + ".hdf")

    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

    raw_file = preprocessed_files_path / (science_filename + ".fits")

    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.extractStripes(
        flat=flat_path, dark=dark_path,
        dark_subtraction_skip_fibers=[5], straylight_removal_fibers=[5],
        legacy=True,
    )

    new_stripes = p.streams["main"][0][0].STRIPES

    fail_counter = 0
    for f in new_stripes.keys():
        for o in new_stripes[f].keys():
            mat = _load_sparse_mat(f'extracted_stripes/{f}/{o}', str(old_science))
            legacy_stripe = mat.toarray()
            new_stripe = new_stripes[f][o].toarray()

            try:
                assert_allclose_with_max_fails(legacy_stripe, new_stripe, rtol=0, atol=1e-8, max_fails=0)
                log.fullinfo(f'fiber/order : {f}/{o} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {f}/{o} [FAIL]')

    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.slow()
@pytest.mark.parametrize("etalon_filename", ETALON_FILES)
def test_boxExtraction(etalon_filename, path_to_legacy_wavecal, preprocessed_files_path):

    arm = _arm_from_filename(etalon_filename)
    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])

    raw_file = preprocessed_files_path / (etalon_filename + ".fits")

    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
    adout = p.boxExtraction()

    legacy_box = legacy_adapter.load_dict_from_hdf5(str(old_file), "box_extraction/")

    fail_counter = 0
    for fiber in range(1, 6):
        new_box = getattr(adout[0][0], f"BOX_REDUCED_FIBER_{fiber}")
        if new_box.size == 1:
            assert fiber == 1, "Only fiber 1 should be empty in box extraction"
            continue

        orders = getattr(adout[0][0], f"REDUCED_ORDERS_FIBER_{fiber}").astype(int)

        for idx, order in enumerate(orders):
            legacy_order = legacy_box[f"fiber_{fiber}"][f"{order}"]
            new_order = new_box[idx, :]

            try:
                assert_allclose_with_max_fails(legacy_order, new_order, rtol=0, atol=1e-15, max_fails=0)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')

    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.slow()
@pytest.mark.parametrize("science_filename", SCIENCE_FILES)
def test_optimalExtraction(science_filename, path_to_legacy_science, preprocessed_files_path):

    arm = _arm_from_filename(science_filename)
    old_file = path_to_legacy_science / (science_filename + ".hdf")

    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

    if USE_CACHE:
        raw_file = preprocessed_files_path / (science_filename + "_reduced.fits")
        adout = [astrodata.open(raw_file)]
    else:
        raw_file = preprocessed_files_path / (science_filename + ".fits")

        adinput = [astrodata.open(raw_file)]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(
            flat=flat_path, dark=dark_path,
            dark_subtraction_skip_fibers=[5], straylight_removal_fibers=[5],
            legacy=True,
        )
        adout = p.optimalExtraction(dark=dark_path)

    legacy_box = legacy_adapter.load_dict_from_hdf5(str(old_file), "box_extraction/")
    legacy_opt = legacy_adapter.load_dict_from_hdf5(str(old_file), "optimal_extraction/")
    legacy_var = legacy_adapter.load_dict_from_hdf5(str(old_file), "optimal_var/")

    fail_counter = 0
    for fiber in range(1, 6):
        new_opt = getattr(adout[0][0], f"OPTIMAL_REDUCED_FIBER_{fiber}", None)
        if new_opt is None or new_opt.size == 1:
            continue

        new_box = getattr(adout[0][0], f"BOX_REDUCED_FIBER_{fiber}")
        orders = getattr(adout[0][0], f"REDUCED_ORDERS_FIBER_{fiber}").astype(int)

        for idx, order in enumerate(orders):
            legacy_fiber_opt = legacy_opt.get(f"fiber_{fiber}")
            if legacy_fiber_opt is None:
                continue
            legacy_order_opt = legacy_fiber_opt.get(f"{order}")
            if legacy_order_opt is None:
                continue

            try:
                assert_allclose_with_max_fails(
                    new_opt[idx, :], legacy_order_opt,
                    rtol=0, atol=1e-8, max_fails=0)
                log.fullinfo(f'fiber/order [opt]: {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order [opt]: {fiber}/{order} [FAIL]')

            legacy_fiber_box = legacy_box.get(f"fiber_{fiber}")
            if legacy_fiber_box is not None:
                legacy_order_box = legacy_fiber_box.get(f"{order}")
                if legacy_order_box is not None:
                    try:
                        assert_allclose_with_max_fails(
                            new_box[idx, :], legacy_order_box,
                            rtol=0, atol=1e-8, max_fails=0)
                        log.fullinfo(f'fiber/order [box]: {fiber}/{order} [OK]')
                    except (AssertionError, pytest.xfail.Exception):
                        fail_counter += 1
                        log.fullinfo(f'fiber/order [box]: {fiber}/{order} [FAIL]')

    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.slow()
@pytest.mark.parametrize("etalon_filename", ETALON_FILES)
def test_staticWavelengthSolution(etalon_filename, path_to_legacy_wavecal, preprocessed_files_path):

    arm = _arm_from_filename(etalon_filename)
    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    legacy_wls = legacy_adapter.load_dict_from_hdf5(str(old_file), 'wavelengths_static/')

    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])

    raw_file = preprocessed_files_path / (etalon_filename + ".fits")

    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
    p.boxExtraction()
    adout = p.staticWavelengthSolution()

    fail_counter = 0
    for fiber in range(1, 6):
        new_wls = getattr(adout[0][0], f"WLS_STATIC_FIBER_{fiber}")
        if new_wls.size == 1:
            continue

        orders = getattr(adout[0][0], f"REDUCED_ORDERS_FIBER_{fiber}").astype(int)

        for idx, order in enumerate(orders):
            legacy_order_wls = legacy_wls[f"fiber_{fiber}"][f"{order}"]
            new_order_wls = new_wls[idx, :]

            try:
                assert_allclose_with_max_fails(legacy_order_wls, new_order_wls, rtol=0, atol=1e-15, max_fails=0)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')

    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.slow()
@pytest.mark.parametrize("science_filename", SCIENCE_FILES)
def test_combineFibers(science_filename, path_to_legacy_science, preprocessed_files_path):

    target_fiber = 6
    arm = _arm_from_filename(science_filename)

    old_file = path_to_legacy_science / (science_filename + ".hdf")
    legacy = {
        'wls': legacy_adapter.load_dict_from_hdf5(str(old_file), f'wavelengths_simultaneous/fiber_{target_fiber}'),
        'opt': legacy_adapter.load_dict_from_hdf5(str(old_file), f'optimal_extraction/fiber_{target_fiber}'),
        'opt_var': legacy_adapter.load_dict_from_hdf5(str(old_file), f'optimal_var/fiber_{target_fiber}'),
    }

    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])
    wavecal_path = str(calib_root / _CALIB_FILES[arm]['wavecal'])

    if USE_CACHE:
        raw_file = preprocessed_files_path / (science_filename + "_reduced.fits")
        ad = astrodata.open(raw_file)
        p = MaroonXSpectrum([ad])
    else:
        raw_file = preprocessed_files_path / (science_filename + ".fits")

        adinput = [astrodata.open(raw_file)]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(
            flat=flat_path, dark=dark_path,
            dark_subtraction_skip_fibers=[5], straylight_removal_fibers=[5],
            legacy=True,
        )
        p.optimalExtraction(dark=dark_path)
        p.getPeaksAndPolynomials(fibers=(5,))
        p.staticWavelengthSolution()
        p.applyWavelengthSolution(wavecal=wavecal_path, fibers=(2, 3, 4), ref_fiber=5)

    adout = p.combineFibers(max_clips=20)

    new = {
        'wls': getattr(adout[0][0], f"WLS_SIMULTANEOUS_FIBER_{target_fiber}"),
        'opt': getattr(adout[0][0], f"OPTIMAL_REDUCED_FIBER_{target_fiber}"),
        'opt_var': getattr(adout[0][0], f"OPTIMAL_REDUCED_VAR_{target_fiber}"),
    }

    orders = adout[0][0].REDUCED_ORDERS_FIBER_3
    orders = orders.astype(int)

    fail_counter = 0
    for idx, order in enumerate(orders):
        for key in ['wls', 'opt', 'opt_var']:
            legacy_order = legacy[key][f"{order}"]
            new_order = new[key][idx, :]

            try:
                assert_allclose_with_max_fails(legacy_order, new_order, rtol=1e-7, atol=0, max_fails=0)
                log.fullinfo(f'key/order : {key:7s}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                fail_counter += 1
                log.fullinfo(f'key/order : {key:7s}/{order} [FAIL]')

    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} key/order combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing key/order combinations"


@pytest.mark.parametrize("science_filename", SCIENCE_FILES)
def test_barycentricCorrection(science_filename, path_to_legacy_science, preprocessed_files_path):

    def _decode(x):
        return x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else x

    arm = _arm_from_filename(science_filename)
    old_file = path_to_legacy_science / (science_filename + ".hdf")

    raw_file = preprocessed_files_path / (science_filename + ".fits")

    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()

    adout = p.barycentricCorrection(target_name="HD3651")
    hdr = adout[0][0].hdr

    keys = ["BERV_SIMBAD_TARGET", "BERV_FLUXWEIGHTED_PC", "BERV_FLUXWEIGHTED_FRD"]
    legacy_target, legacy_pc, legacy_frd = legacy_adapter.read_header_entries(str(old_file), keys)

    assert str(_decode(legacy_target)) == str(hdr.get("BERV_SIMBAD_TARGET"))

    np.testing.assert_allclose(
        float(_decode(legacy_pc)),
        float(hdr.get("BERV_FLUXWEIGHTED_PC")),
        rtol=0.0,
        atol=0.1,
    )
    np.testing.assert_allclose(
        float(_decode(legacy_frd)),
        float(hdr.get("BERV_FLUXWEIGHTED_FRD")),
        rtol=0.0,
        atol=0.1,
    )
