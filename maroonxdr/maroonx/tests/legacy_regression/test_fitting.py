import astrodata
import numpy as np
import pandas as pd
import pytest
from gempy.utils import logutils

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

from . import legacy_adapter

_CALIB_FILES = {
    'RED': {
        'flat': 'processed_flat/20250701T171955Z_DDDDF_r_0002_DFFFF_flat.fits',
    },
    'BLUE': {
        'flat': 'processed_flat/20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits',
    },
}

logutils.config(file_name="test_fitting_v2.log", mode="debug", stomp=True)
log = logutils.get_logger("test_fitting_v2")
log.setLevel("DEBUG")

# =========================================================
# TESTS
# =========================================================
USE_CACHE = True

ETALON_FILES = [
    '20250717T163124Z_DEEEE_b_0010',
    '20250717T163124Z_DEEEE_r_0004',
]


def _arm_from_filename(filename):
    return 'BLUE' if '_b_' in filename else 'RED'


@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize("etalon_filename", ETALON_FILES)
def test_getPeaksAndPolynomials(path_to_legacy_wavecal, preprocessed_files_path, etalon_filename):

    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    assert "fiber" in old_peak_data.columns

    if USE_CACHE:
        adout = [astrodata.open(str(preprocessed_files_path / (etalon_filename + "_wavecal.fits")))]
    else:
        arm = _arm_from_filename(etalon_filename)
        calib_root = preprocessed_files_path / 'calibrations'
        flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])

        adinput = [astrodata.open(str(preprocessed_files_path / (etalon_filename + ".fits")))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
        p.boxExtraction()
        adout = p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))

    new_peak_data = adout[0][0].PEAKS.to_pandas()

    new_peak_data.sort_values(by=['FIBER', 'ORDER', 'CENTER'], inplace=True)
    old_peak_data.sort_values(by=['fiber', 'order', 'center'], inplace=True)

    assert old_peak_data.shape == new_peak_data.shape

    assert new_peak_data["FIBER"].value_counts().equals(old_peak_data["fiber"].value_counts())

    col_tolerances = {
        'amplitude': {'atol': 1e-5, 'rtol': 0},
        'center':    {'atol': 1e-5, 'rtol': 0},
        'fiber':     {'atol': 0,    'rtol': 0},
        'order':     {'atol': 0,    'rtol': 0},
        'offset':    {'atol': 1e-5, 'rtol': 0},
        'sigma1':    {'atol': 1e-5, 'rtol': 0},
        'sigma2':    {'atol': 1e-5, 'rtol': 0},
        'width':     {'atol': 1e-5, 'rtol': 0},
    }

    for col, tol in col_tolerances.items():
        try:
            np.testing.assert_allclose(
                new_peak_data[col.upper()].values,
                old_peak_data[col.lower()].values,
                **tol,
                err_msg=f"Column '{col}' mismatch",
            )
        except AssertionError as err:
            log.fullinfo(f"Column {col} mismatch: {err}")
            raise err


@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize("etalon_filename", ETALON_FILES)
def test_fitAndApplyEtalonWls(path_to_legacy_wavecal, preprocessed_files_path, etalon_filename):

    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    old_peak_data = pd.read_hdf(old_file, '/peak_data').reset_index()
    assert "fiber" in old_peak_data.columns

    if USE_CACHE:
        adout = [astrodata.open(str(preprocessed_files_path / (etalon_filename + "_wavecal.fits")))]
    else:
        arm = _arm_from_filename(etalon_filename)
        calib_root = preprocessed_files_path / 'calibrations'
        flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])

        adinput = [astrodata.open(str(preprocessed_files_path / (etalon_filename + ".fits")))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
        p.boxExtraction()
        p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))
        p.staticWavelengthSolution()
        adout = p.fitAndApplyEtalonWls()

    new_peak_data = adout[0][0].PEAK_DATA.to_pandas()

    new_cols = {c.lower() for c in new_peak_data.columns}
    old_cols = {c.lower() for c in old_peak_data.columns}
    assert new_cols.issubset(old_cols)

    assert new_peak_data["FIBER"].value_counts().equals(old_peak_data["fiber"].value_counts())

    col_tolerances = {
        'amplitude':          {'atol': 1e-5, 'rtol': 0},
        'center':             {'atol': 1e-5, 'rtol': 0},
        'fiber':              {'atol': 0,    'rtol': 0},
        'order':              {'atol': 0,    'rtol': 0},
        'offset':             {'atol': 1e-5, 'rtol': 0},
        'sigma1':             {'atol': 1e-5, 'rtol': 0},
        'sigma2':             {'atol': 1e-5, 'rtol': 0},
        'width':              {'atol': 1e-5, 'rtol': 0},
        'wavelength_by_thar': {'atol': 1e-8, 'rtol': 0},
        'm':                  {'atol': 0,    'rtol': 0},
        'm_fraction':         {'atol': 1e-5, 'rtol': 0},
        'dispersion_mps':     {'atol': 1e-5, 'rtol': 0},
    }

    for col, tol in col_tolerances.items():
        if col.upper() not in new_peak_data.columns:
            continue
        if col.lower() not in old_peak_data.columns:
            continue
        try:
            np.testing.assert_allclose(
                new_peak_data[col.upper()].values,
                old_peak_data[col.lower()].values,
                **tol,
                err_msg=f"Column '{col}' mismatch",
            )
        except AssertionError as err:
            log.fullinfo(f"Column {col} mismatch: {err}")
            raise err


@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize("etalon_filename", ETALON_FILES)
def test_dynamicWavelengthSolution(path_to_legacy_wavecal, preprocessed_files_path, etalon_filename):

    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    legacy_wls = legacy_adapter.load_dict_from_hdf5(str(old_file), 'wavelengths_dynamic/')

    if USE_CACHE:
        adout = [astrodata.open(str(preprocessed_files_path / (etalon_filename + "_wavecal.fits")))]
    else:
        arm = _arm_from_filename(etalon_filename)
        calib_root = preprocessed_files_path / 'calibrations'
        flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])

        adinput = [astrodata.open(str(preprocessed_files_path / (etalon_filename + ".fits")))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
        p.boxExtraction()
        p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))
        p.staticWavelengthSolution()
        adout = p.fitAndApplyEtalonWls()

    fail_counter = 0
    for fiber in range(1, 6):
        new_wls = getattr(adout[0][0], f"WLS_DYNAMIC_FIBER_{fiber}")
        if new_wls.size == 1:
            continue

        orders = getattr(adout[0][0], f"REDUCED_ORDERS_FIBER_{fiber}").astype(int)

        for idx, order in enumerate(orders):
            legacy_order_wls = legacy_wls[f"fiber_{fiber}"][f"{order}"]
            new_order_wls = new_wls[idx, :]

            try:
                np.testing.assert_allclose(legacy_order_wls, new_order_wls, rtol=1e-8)
                log.fullinfo(f'fiber/order : {fiber}/{order} [OK]')
            except AssertionError:
                fail_counter += 1
                log.fullinfo(f'fiber/order : {fiber}/{order} [FAIL]')
    assert fail_counter == 0, f"{fail_counter} failing f/o combinations"


@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize("etalon_filename", ETALON_FILES)
def test_fiber_drifts_ETALON(path_to_legacy_wavecal, preprocessed_files_path, etalon_filename):

    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    hdr_keys = [
        'Drift_Fiber1', 'Drift_Fiber2',
        'Drift_Fiber3', 'Drift_Fiber4', 'Instrument_Drift',
    ]
    old_header_drifts = legacy_adapter.read_header_entries(str(old_file), hdr_keys)

    legacy_drifts = {}
    for fiber, raw_entry in enumerate(old_header_drifts, start=1):
        drift = float(raw_entry.decode().split()[0]) if raw_entry is not None else None
        legacy_drifts[f"fiber_{fiber}"] = drift

    if USE_CACHE:
        adout = [astrodata.open(str(preprocessed_files_path / (etalon_filename + "_wavecal.fits")))]
    else:
        arm = _arm_from_filename(etalon_filename)
        calib_root = preprocessed_files_path / 'calibrations'
        flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])

        adinput = [astrodata.open(str(preprocessed_files_path / (etalon_filename + ".fits")))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
        p.boxExtraction()
        p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))
        p.staticWavelengthSolution()
        adout = p.fitAndApplyEtalonWls()

    drifts = adout[0].fiber_drifts()

    precission = 0.1  # m/s

    for fiber, drift in enumerate(drifts, start=1):
        legacy_drift = legacy_drifts[f"fiber_{fiber}"]
        if legacy_drift is None:
            assert drift is None, f"Fiber {fiber} drift mismatch: legacy is None, but new is {drift}"
            continue

        assert abs(legacy_drift - drift) < precission, f"Fiber {fiber} drift mismatch: {legacy_drift:.2f} vs {drift:.2f} [m/s]"
