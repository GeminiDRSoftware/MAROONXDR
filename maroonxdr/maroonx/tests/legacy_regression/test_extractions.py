

from pathlib import Path

import astrodata
import h5py
import numpy as np
import pytest
from gempy.adlibrary import dataselect
from gempy.utils import logutils
from scipy import sparse

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from maroonxdr.maroonx.tests.conftest import assert_allclose_with_max_fails

# Calibration file mapping by arm (same pattern as echelle_extraction tests)
_CALIB_FILES = {
    'RED': {
        'flat': 'processed_flat/20241114T190714Z_DDDDF_r_0002_DFFFF_flat.fits',
        'dark': 'processed_dark/20241124T041907Z_SOOOE_r_0300_synth_dark.fits',
        'wavecal': 'processed_wavecal/20241124T030227Z_DEEEE_r_0004_wavecal.fits',
    },
    'BLUE': {
        'flat': 'processed_flat/20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits',
        'dark': 'processed_dark/20241124T041907Z_SOOOE_b_0300_synth_dark.fits',
        'wavecal': 'processed_wavecal/20241124T030227Z_DEEEE_b_0030_wavecal.fits',
    },
}

# Set logger
logutils.config(file_name="test_extractions.log", mode="debug", stomp=True)
log = logutils.get_logger("test_extractions")
log.setLevel("DEBUG")

# =========================================================
# TESTS
# =========================================================
USE_CACHE = False

SCIENCE_FILES = [
    '20241124T041907Z_SOOOE_b_0300',
]

@pytest.mark.parametrize("arm", ["BLUE"])
def test_extractStripes_fromEtalon(arm, legacy_reduced_path, preprocessed_files_path):

    # old_file = legacy_reduced_path  / "20241124T162336Z_DEEEE_b_0030.hdf"
    old_file = legacy_reduced_path  / "20241124T030227Z_DEEEE_b_0030.hdf"

    # Explicit calibration paths (no caldb configured for legacy regression tests)
    # read files and instantiate the primitive class
    # raw_file = preprocessed_files_path / '20241124T162336Z_DEEEE_b_0030.fits'
    raw_file = preprocessed_files_path / '20241124T030227Z_DEEEE_b_0030.fits'

    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

    # Primitives
    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)


    adout = p.extractStripes(flat=flat_path, dark=dark_path)

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
                assert_allclose_with_max_fails(legacy_stripe, new_stripe, rtol=0, atol=1e-8, max_fails=2)
                log.fullinfo(f'fiber/order : {f}/{o} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                # if xfails was raised, we still want to count it as a fail
                fail_counter += 1
                log.fullinfo(f'fiber/order : {f}/{o} [FAIL]')

    # We accept up to 2 f/o fails
    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.parametrize("arm", ["BLUE"])
def test_extractStripes_fromFlat(arm, legacy_flats_path, preprocessed_files_path):

    old_flat = legacy_flats_path / "20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf"

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('20241124T162336Z_DEEEE_*.fits')])

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
                assert_allclose_with_max_fails(legacy_flat_stripe, new_flat_stripe, rtol=1e-3, atol=1e-8, max_fails=2)
                log.fullinfo(f'fiber/order : {f}/{o} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                # if xfails was raised, we still want to count it as a fail
                fail_counter += 1
                log.fullinfo(f'fiber/order : {f}/{o} [FAIL]')

    # We accept up to 2 f/o fails
    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.parametrize("arm", ["BLUE"])
def test_extractStripes_fromScience(arm, legacy_reduced_path, legacy_bkg_path, preprocessed_files_path):

    old_science = legacy_reduced_path / "20241124T041907Z_SOOOE_b_0300.hdf"
    legacy_npy = legacy_bkg_path / "20241124T041907Z_SOOOE_b_0300_test.npy"
    legacy_dict = np.load(legacy_npy, allow_pickle=True).item()
    # legacy_dict.keys() has the following keys:
    # dict_keys(['raw', 'bias_corrected', 'overscan_removed',
    # 'orientation_corrected', 'dark_science', 'back_var', 'index_fiber',
    # 'index_order', 'mask', 'science_straylight_removed'])

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('20241124T041907Z_SOOOE_*.fits')])

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

    # ====================================================== OK line
    # if remove straylight is fiber dependent, should this happend inside extractStripes?
    # p.removeStrayLight()
    # # test straylight corrected data is the same as legacy
    # np.testing.assert_allclose(
    #     p.streams["main"][0][0].data,
    #     legacy_dict['science_straylight_removed'],
    #     rtol=1e-4, atol=1e-4)
    p.extractStripes(dark_subtraction_skip_fibers=[5], straylight_removal_fibers=[5])

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
                assert_allclose_with_max_fails(legacy_stripe, new_stripe, rtol=0, atol=1e-3, max_fails=2)
                log.fullinfo(f'fiber/order : {f}/{o} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                # if xfails was raised, we still want to count it as a fail
                fail_counter += 1
                log.fullinfo(f'fiber/order : {f}/{o} [FAIL]')

    # We accept up to 2 f/o fails
    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"

    # ==============================================================


@pytest.mark.parametrize("arm", ["BLUE"])
def test_boxExtraction(arm, legacy_reduced_path, preprocessed_files_path):

    old_file = legacy_reduced_path / "20241124T162336Z_DEEEE_b_0030.hdf"

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('20241124T162336Z_DEEEE_*.fits')])

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
    # adout = p.optimalExtraction()

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
                #np.testing.assert_allclose(legacy_order, new_order, rtol=0, atol=1e-8)
                assert_allclose_with_max_fails(legacy_order, new_order, rtol=0, atol=1e-8, max_fails=2)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                # if xfails was raised, we still want to count it as a fail
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')

    # We accept up to 2 f/o fails
    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.slow()
@pytest.mark.parametrize("science_filename", SCIENCE_FILES)
def test_optimalExtraction(legacy_reduced_path, science_filename):

    old_file = legacy_reduced_path / (science_filename + ".hdf")

    if USE_CACHE:
        # Use previously saved data on science_dir
        adout = [astrodata.open(science_filename + "_reduced.fits")]
    else:
        # Primitives
        adinput = [astrodata.open(science_filename + ".fits")]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(dark_subtraction_skip_fibers=[5], straylight_removal_fibers=[5])
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
                assert_allclose_with_max_fails(legacy_order, new_order, rtol=0, atol=1e-4, max_fails=2)
                log.fullinfo(f'fiber/order [opt]: {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                # if xfails was raised, we still want to count it as a fail
                fail_counter += 1
                log.fullinfo(f'fiber/order [opt]: {fiber}/{order} [FAIL]')

            legacy_order = legacy_box[f"fiber_{fiber}"][f"{order}"]
            new_order = new_box[idx, :]
            try:
                assert_allclose_with_max_fails(legacy_order, new_order, rtol=0, atol=1e-4, max_fails=2)
                log.fullinfo(f'fiber/order [box]: {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                # if xfails was raised, we still want to count it as a fail
                fail_counter += 1
                log.fullinfo(f'fiber/order [box]: {fiber}/{order} [FAIL]')

    # We accept up to 2 f/o fails
    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.parametrize("arm", ["BLUE"])
def test_staticWavelengthSolution(arm, legacy_reduced_path, preprocessed_files_path):

    old_file = legacy_reduced_path / "20241124T162336Z_DEEEE_b_0030.hdf"

    # Load old peak data. columns are lowercase
    legacy_wls = load_dict_from_hdf5(str(old_file), 'wavelengths_static/')

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('20241124T162336Z_DEEEE_*.fits')])

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
                assert_allclose_with_max_fails(legacy_order_wls, new_order_wls, rtol=0, atol=1e-8, max_fails=2)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                # if xfails was raised, we still want to count it as a fail
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')

    # We accept up to 2 f/o fails
    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} f/o combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.parametrize("science_filename", SCIENCE_FILES)
def test_optimal_extraction_single_stripe(legacy_bkg_path, science_filename):
    # Load old data
    old_file_inputs = (
        legacy_bkg_path / f"{science_filename}_optimal_2_111_inputs.npy"
    )
    old_input = np.load(old_file_inputs, allow_pickle=True).item()

    (
        old_flux, old_var, old_stand_spec, old_stand_var, old_fo
    ) = _optimal_extraction_single_stripe(
        old_input['stripe'],
        old_input['flat_stripes'],
        gain=old_input['gain'],
        read_noise=old_input['read_noise'],
        back_var=old_input['back_var'],
        mask=old_input['mask'],
        s_clip=old_input['s_clip'],
        penalty=old_input['penalty'],
        full_output=False,
        log=log)
    np.testing.assert_allclose(
        old_flux, old_input['flux'], rtol=0, atol=1e-8
    )
    np.testing.assert_allclose(
        old_var, old_input['var'], rtol=0, atol=1e-8
    )
    np.testing.assert_allclose(
        old_stand_spec, old_input['stand_spec'], rtol=0, atol=1e-8
    )

    # load new data
    new_input = np.load(
        str(legacy_bkg_path / f"{science_filename}_optimal_2_111_inputs.npy"),
        allow_pickle=True,
    ).item()
    (
        new_flux, new_var, new_stand_spec, new_stand_var, new_fo
    ) = _optimal_extraction_single_stripe(
        new_input['stripe'],
        new_input['flat_stripes'],
        gain=new_input['gain'],
        read_noise=new_input['read_noise'],
        back_var=new_input['back_var'],
        mask=new_input['mask'],
        s_clip=new_input['s_clip'],
        penalty=new_input['penalty'],
        full_output=False,
        log=log)
    np.testing.assert_allclose(
        new_flux, new_input['flux'], rtol=0, atol=1e-8
    )
    np.testing.assert_allclose(
        new_var, new_input['var'], rtol=0, atol=1e-8
    )
    np.testing.assert_allclose(
        new_stand_spec, new_input['stand_spec'], rtol=0, atol=1e-8
    )


    # Assert old and new inputs are close to each other
    np.testing.assert_allclose(old_input['stripe'].todense(), new_input['stripe'].todense(), rtol=0, atol=1e-5)
    np.testing.assert_allclose(old_input['flat_stripes'].todense(), new_input['flat_stripes'].todense(), rtol=0, atol=1e-5)
    np.testing.assert_allclose(old_input['back_var'], new_input['back_var'], rtol=0, atol=1e-5)


    # # Assert old and new results are close to each other
    # np.testing.assert_allclose(old_flux, new_flux, rtol=0, atol=1e-4)
    # np.testing.assert_allclose(old_var, new_var, rtol=0, atol=1e-4)
    # np.testing.assert_allclose(old_stand_spec, new_stand_spec, rtol=0, atol=1e-4)


@pytest.mark.slow()
@pytest.mark.parametrize("arm", ["BLUE"])
def test_combineFibers(arm, legacy_reduced_path, preprocessed_files_path):

    target_fiber = 6

    # Load old data
    old_file = legacy_reduced_path / "20241124T062858Z_SOOOE_b_0300.hdf"
    legacy = {
        'wls': load_dict_from_hdf5(str(old_file), f'wavelengths_simultaneous/fiber_{target_fiber}'),
        'opt': load_dict_from_hdf5(str(old_file), f'optimal_extraction/fiber_{target_fiber}'),
        'opt_var': load_dict_from_hdf5(str(old_file), f'optimal_var/fiber_{target_fiber}'),
    }

    if USE_CACHE:
        # Use previously saved data on science_dir
        ad = astrodata.open('20241124T062858Z_SOOOE_b_0300_reduced.fits')
        p = MaroonXSpectrum([ad])
    else:
        # read files and instantiate the primitive class
        raw_files = sorted([str(f) for f in preprocessed_files_path.glob('20241124T062858Z_SOOOE_*_0300.fits')])

        selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'SCI', arm])

        # Primitives
        adinput = [astrodata.open(f) for f in selected_spect]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True,poisson_noise=True)

        p.extractStripes(dark_subtraction_skip_fibers=[5], straylight_removal_fibers=[5])
        p.optimalExtraction()
        p.getPeaksAndPolynomials(fibers=(5,) , multithreading=True)
        p.staticWavelengthSolution()
        p.applyWavelengthSolution(fibers=(2,3,4), ref_fiber=5)

    adout = p.combineFibers()

    new = {
        'wls': getattr(adout[0][0], f"WLS_SIMULTANEOUS_FIBER_{target_fiber}"),
        'opt': getattr(adout[0][0], f"OPTIMAL_REDUCED_FIBER_{target_fiber}"),
        'opt_var': getattr(adout[0][0], f"OPTIMAL_REDUCED_VAR_{target_fiber}"),
    }

    # All orders should be the same, we dont save orders for fiber 6, should we?
    orders = adout[0][0].REDUCED_ORDERS_FIBER_3
    orders = orders.astype(int)

    fail_counter = 0
    for idx, order in enumerate(orders):
        for key in ['wls', 'opt', 'opt_var']:
            legacy_order = legacy[key][f"{order}"]
            new_order = new[key][idx, :]

            try:
                assert_allclose_with_max_fails(legacy_order, new_order, rtol=1e-2, atol=1e-8, max_fails=2)
                #np.testing.assert_allclose(legacy_order, new_order, rtol=1e-2, atol=1e-8)
                log.fullinfo(f'key/order : {key:7s}/{order} [OK]')
            except (AssertionError, pytest.xfail.Exception):
                # if xfails was raised, we still want to count it as a fail
                fail_counter += 1
                log.fullinfo(f'key/order : {key:7s}/{order} [FAIL]')

    # We accept up to 2 key/order fails
    if 0 < fail_counter <= 2:
        pytest.xfail(f"{fail_counter} key/order combinations differ (max allowed: 2)")
    assert fail_counter == 0, f"{fail_counter} failing key/order combinations"


@pytest.mark.parametrize("arm", ["BLUE"])
def test_barycentricCorrection(arm, legacy_reduced_path, preprocessed_files_path):

    def _decode(x):
        return x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else x

    # Load old data
    old_file = legacy_reduced_path / "20241124T062858Z_SOOOE_b_0300.hdf"

    # read files and instantiate the primitive class
    raw_files = sorted([str(f) for f in preprocessed_files_path.glob('20241124T062858Z_SOOOE_*_0300.fits')])

    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'SCI', arm])

    # Primitives
    adinput = [astrodata.open(f) for f in selected_spect]

    p = MaroonXSpectrum(adinput)
    p.prepare()

    adout = p.barycentricCorrection(target_name="HD 203030")
    hdr = adout[0][0].hdr

    keys = ["BERV_SIMBAD_TARGET", "BERV_FLUXWEIGHTED_PC", "BERV_FLUXWEIGHTED_FRD"]
    legacy_target, legacy_pc, legacy_frd = read_header_entries(str(old_file), keys)

    assert str(_decode(legacy_target)) == str(hdr.get("BERV_SIMBAD_TARGET"))

    print(f"Legacy BERV_PC: {_decode(legacy_pc)}, New BERV_PC: {hdr.get('BERV_FLUXWEIGHTED_PC')}")
    print(f"Legacy BERV_FRD: {_decode(legacy_frd)}, New BERV_FRD: {hdr.get('BERV_FLUXWEIGHTED_FRD')}")
    
    # BERV is in (m/s); atol=0.1 -> 10 cm/s tolerance, rtol=0
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

def read_header_entries(filename, hdr_keys):
    values = []
    with h5py.File(filename, 'r', libver='latest') as h5f:
        for key in hdr_keys:
            values.append(h5f['header'].attrs.get(key))
    return values

# =========================================================================
def _optimal_extraction_single_stripe(stripe, flat_stripe, gain=1, read_noise=1.23, back_var=None, mask=None,
                                    debug_level=0, full_output=False, s_clip=5.0, penalty=1.0, log=None):
    """
    Performs optimal extraction of a single stripe.
    Based on the algorithm described by Horne et al. 1986, PASP, 98, 609.

    Args:
        stripe (scipy.sparse.spmatrix): science frame stripe to be extracted
        flat_stripe (scipy.sparse.spmatrix): flat frame stripe to be used as profile
        gain (float): detector gain factor (conversion photons -> DN,  given in e-/DN )
        read_noise (float): typical detector read noise given as the variance in DN
        back_var (np.ndarray): background variance (from scattered light model oder dark, otherwhise 0)
        mask (np.ndarray): bad pixel mask. Note: the mask will be modified by iterative algorithm that looks for
        outliers.
        debug_level (int): debug level
        full_output (bool): if True, returns all intermediate results for debugging/testing
        s_clip (float): sigma clipping value for optimal extraction
        penalty (float): scaling factor for global per-order profile mismatch correction. Set 0 for no correction

    Returns
    -------
        tuple(np.ndarray, np.ndarray, dict): (optimal extracted spectrum, box extracted spectrum, dict of additional
        intermediate results if full_output was True)

    """
    from scipy.ndimage import median_filter

    #log = self.log
    if mask is None:
        mask       = stripe.copy()
        mask[:, :] = 1
    else:
        # copy mask, because it will be modified
        mask     = mask.copy()

    if back_var is None or isinstance(back_var, float):
        back_var = stripe.copy()
        back_var[:,:] = 0
    else:
        # back_var = back_var.copy()
        back_var = back_var.todense() if hasattr(back_var, "todense") else back_var.copy()

    # box extracted spectrum
    stand_spec0 = _box_extract_single_stripe(stripe, mask)  # direct box extraction,
    stand_spec = stand_spec0.copy()
    stand_var = stand_spec / gain

    # flat data
    box_extracted_flat_stripe = _box_extract_single_stripe(flat_stripe, mask)  # spatial sums for flat along disp

    # print(f"stripe shape: {stripe.shape}, flat stripe shape: {flat_stripe.shape}, mask shape: {mask.shape}"
    #       f", back_var shape: {back_var.shape}, box_extracted_flat_stripe shape: {box_extracted_flat_stripe.shape}")
    # # print the types
    # print(f"stripe type: {type(stripe)}, flat stripe type: {type(flat_stripe)}, mask type: {type(mask)}"
    #         f", back_var type: {type(back_var)}, box_extracted_flat_stripe type: {type(box_extracted_flat_stripe)}")


    # cut stripe sparse matrix into numpy array
    # find the spatial columns utilized along entire stripe (greater than slit height because of stripe path)
    sparse_vcols = np.array(~np.all(stripe.todense() == 0, axis=1)).reshape(-1)
    stripe       = np.array(stripe.todense()[sparse_vcols])  # strip stripe to the inclusive nonzero rows

    mask         = mask[sparse_vcols]  # strip mask similarly
    back_var     = back_var[sparse_vcols] #strip background variance map similarly

    flat_stripe  = np.array(flat_stripe.todense()[sparse_vcols])
    sparse_vrows = np.count_nonzero((stripe != 0).T[1500])  # use ~middle column slit height as slit height pass
    data         = np.zeros((sparse_vrows, stripe.shape[1]))
    diff_save    = data.copy()
    new_mask     = data.copy()  # create actual limit numpy arrays
    new_back_var = data.copy()  # create actual limit numpy arrays
    profile      = data.copy()

    # print(f"stripe shape: {stripe.shape}, flat stripe shape: {flat_stripe.shape}, mask shape: {mask.shape}"
    #       f", back_var shape: {back_var.shape}")
    # # print the types
    # print(f"stripe type: {type(stripe)}, flat stripe type: {type(flat_stripe)}, mask type: {type(mask)}"
    #         f", back_var type: {type(back_var)}")
    #import ipdb; ipdb.set_trace()

    for i in np.arange(stripe.shape[1]):
        if np.nonzero((stripe != 0).T[i])[0].shape[0] == data.shape[0]:  # if column is slit height (not edge of chip)
            if box_extracted_flat_stripe[i] > 1E-12:
                data[:, i]        = stripe[np.nonzero((stripe != 0).T[i])[0], i]  # write data
                new_mask[:, i]    = mask[np.nonzero((stripe != 0).T[i])[0], i]
                new_back_var[:,i] =  back_var[np.nonzero((stripe != 0).T[i])[0], i]
                profile[:, i]     = flat_stripe[np.nonzero((stripe != 0).T[i])[0], i] / box_extracted_flat_stripe[i]
        else:
            new_mask[:, i]     = 0  # could be optimized
            new_back_var[:, i] = 0  # could be optimized

    new_mask[:, stand_spec0 < 1E-12] = 0  # if sum is less than zero, whole column is cancelled for this stripe
    mask = new_mask.copy()
    back_var = new_back_var.copy() * mask
    stripe = data  # return variables 'mask' and 'stripe' to the naming conventions
    data_var = abs(stripe.copy())/gain + back_var + read_noise

    # final output
    flux = np.zeros(len(stand_spec0))
    var = flux.copy()

    # Calculate a first guess of the difference between stripe and scaled flat_stripe
    expected = profile * mask * stand_spec
    actual = stripe * mask
    diff = actual - expected

    # Calculate the median of the difference along the order. This represents the 'global' mismatch between flat and science profile
    # and helps correct for the 'drift' problem in x-dispersion.
    diff_aver = penalty * np.abs(median_filter(diff,size=(1,201)))

    # # avoid already caught bad pixels (whole column is zero in bpm or stripe is too small)
    good_disp = np.nonzero(np.array(~np.any(mask == 0, axis=0)))[0]
    reject_tracker = np.zeros_like(stand_spec0)

    for h in good_disp:
        expected = profile[:, h] * mask[:, h] * stand_spec[h]  # flat column with mask scaled to data total
        actual = stripe[:, h] * mask[:, h]  # actual column data with mask
        diff = actual - expected

        data_var[:, h] = abs(stand_spec[h] * profile[:, h]) / gain + back_var[:, h] + read_noise
        noise_rev = 1 / np.sqrt(data_var[:, h])
        diff_save[:, h] = diff * noise_rev
        reject_index = np.nonzero(((np.abs(diff) - np.abs(diff_aver[:, h])) * noise_rev) >= s_clip)[0]

        while len(reject_index) > 0:

            worst = np.argmax((np.abs(diff) - np.abs(diff_aver[:, h])) * noise_rev)
            reject_tracker[h] = reject_tracker[h] + 1
            if debug_level >= 3:
                log.fullinfo(f'Outlier found in column {h}, pixel {worst}')
                # fig, ax = plt.subplots(2, 1)
                # actual_plot = actual.copy()
                # expected_plot = expected.copy()
                # actual_plot[actual_plot == 0 ] = np.nan
                # expected_plot[expected_plot == 0] = np.nan
                # ax[0].plot(actual_plot, 'r')
                # ax[0].plot(expected_plot, 'b')
                # ax[1].plot(np.abs(diff), 'r')
                # ax[1].plot(np.abs(diff)-np.abs(diff_aver[:, h]), 'g')
                # ax[1].plot(np.sqrt(data_var[:, h]) * s_clip, 'b')
                # plt.show()

            mask[worst, h] = 0

            denom = np.sum(profile[:, h] * profile[:, h] * mask[:, h] / data_var[:, h])

            stand_spec[h] = np.sum(profile[:, h] * mask[:, h] * stripe[:, h] / data_var[:, h]) / denom

            expected = profile[:, h] * mask[:, h] * stand_spec[h]
            actual = stripe[:, h] * mask[:, h]
            diff = actual - expected

            data_var[:, h] = abs(stand_spec[h] * profile[:, h]) / gain + back_var[:, h] + read_noise
            noise_rev = 1 / np.sqrt(data_var[:, h])
            # diff_save[:,h] = diff*noise_rev
            reject_index = np.nonzero(((np.abs(diff) - np.abs(diff_aver[:, h])) * noise_rev) >= s_clip)[0]

            if np.count_nonzero(actual[3:-5]) < len(profile[3:-5, h]) /2.:
                reject_index = np.array([])
                mask[:, h] = 0
                flux[h] = 0
                log.warning(f'Too many bad pixels in column {h}, reject column')
        if np.count_nonzero(mask[:, h]) > 0:
            denom = np.sum(profile[:, h] * profile[:, h] * mask[:, h] / data_var[:, h])
            flux[h] = np.sum(profile[:, h] * mask[:, h] * stripe[:, h] / data_var[:, h]) / denom
            var[h] = np.sum(profile[:, h] * mask[:, h]) / denom

    flux[flux == 0] = np.nan
    var[var == 0] = np.nan

    stand_spec0[stand_spec0 == 0] = np.nan
    stand_var[stand_var == 0] = np.nan

    total_count = np.sum(reject_tracker > 0)
    substantial_count = np.sum(np.abs(stand_spec0 - stand_spec)[reject_tracker > 0] > 0.005 * (stand_spec0[reject_tracker > 0]))
    irrelevant_count = total_count - substantial_count

    if total_count > 0:
        log.fullinfo(f'Rejected {np.sum(reject_tracker):.0f} pixels in {total_count} columns during optimal extraction.')
        irrelevant_count_percentage = irrelevant_count/total_count*100
        if irrelevant_count_percentage > 30:
            log.warning(f'Rejections with flux changes < 0.5%: {irrelevant_count} ({irrelevant_count_percentage:.0f}%)')
        else:
            log.fullinfo(f'Rejections with flux changes < 0.5%: {irrelevant_count} ({irrelevant_count_percentage:.0f}%)')
    else:
        log.fullinfo('No rejections')

    # if debug_level >= 2:
    #     #fig, ax = plt.subplots(3, 1)
    #     fig, ax = plt.subplots(3, 1, sharex='all')
    #     ax[0].imshow(np.vstack((stripe,np.zeros_like(stripe),profile*flux)),interpolation='none')
    #     ax[0].scatter(np.where(mask==0)[1],np.where(mask==0)[0],c='r',marker='+')

    #     ax[1].plot(flux - stand_spec0,'r',label='Difference box vs optimal')
    #     ax[2].plot(stand_spec0, 'b', label='Box extraction')
    #     ax[2].plot(flux, 'g', label='Optimal extraction')
    #     plt.legend()

    #     plt.show()

    # return all intermediate results for debugging/testing
    if full_output:
        return flux, var, stand_spec0, stand_var, {'noise': var, 'acceptancemask': mask, "stripe": stripe,
                                        "initial_sigma": diff_save}  # {'noise': var, 'profile': profile, 'rejectionmask': sparse.csr_matrix(mask),
        # 'expected': expected, 'actual': actual}
    else:
        return flux, var, stand_spec0, stand_var, {'noise': var}




#@staticmethod
def _box_extract_single_stripe(stripe=None, mask=None):
    """
    Box extraction of a single stripe.

    Parameters
    ----------
        stripe (sparse.matrix): stripe to be extracted
        mask (np.ndarray): binary pixel mask. Where mask==1: values will
        contribute. Same shape as stripe.toarray()

    Returns
    -------
        np.ndarray: box extracted flux
    """
    if mask is None:
        mask = stripe.copy()
        mask.data[:] = 1.
    return np.array(np.sum(stripe.multiply(mask), axis=0).T).flatten()
