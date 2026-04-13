

from pathlib import Path

import astrodata
import h5py
import numpy as np
import pandas as pd
import pytest
import tables
from gempy.adlibrary import dataselect
from gempy.utils import logutils

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.maroonx_fit import maroonx_fit, set_logger
from maroonxdr.maroonx.maroonx_utils import load_recordings
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum

from . import legacy_adapter

# Calibration file mapping by arm (same pattern as test_extractions)
_CALIB_FILES = {
    'RED': {
        'flat': 'processed_flat/20241114T190714Z_DDDDF_r_0002_DFFFF_flat.fits',
        'dark': 'processed_dark/20241124T041907Z_SOOOE_r_0300_synth_dark.fits',
    },
    'BLUE': {
        'flat': 'processed_flat/20241114T190714Z_DDDDF_b_0007_DFFFF_flat.fits',
        'dark': 'processed_dark/20241124T041907Z_SOOOE_b_0300_synth_dark.fits',
    },
}

# =========================================================
# FILES TO COMPARE
# =========================================================


# Set logger
logutils.config(file_name="test_fitting.log", mode="debug", stomp=True)
log = logutils.get_logger("test_fitting")
log.setLevel("DEBUG")

# =========================================================
# TESTS
# =========================================================
USE_CACHE = False

ETALONS = [
    # '20241124T030227Z_DEEEE_r_0004',
    '20241124T030227Z_DEEEE_b_0030',
    # '20241124T030040Z_DEEEE_r_0004',
    # '20241124T030436Z_DLLLL_r_0004',
    # '20241124T030623Z_DLLLL_r_0004',
    # '20241124T162149Z_DEEEE_r_0004',
    # '20241124T162336Z_DEEEE_r_0004',
    # '20241124T162540Z_DLLLL_r_0004',
    # '20241124T162727Z_DLLLL_r_0004',
    # '20241124T030040Z_DEEEE_b_0030',
    # '20241124T030436Z_DLLLL_b_0005',
    # '20241124T030623Z_DLLLL_b_0005',
    # '20241124T032549Z_DLLLE_b_0030',
    # '20241124T032715Z_DLLLE_b_0030',
    # '20241124T033308Z_DLLLD_b_0030',
    # '20241124T162149Z_DEEEE_b_0030',
    # '20241124T162336Z_DEEEE_b_0030',
    # '20241124T162540Z_DLLLL_b_0005',
    # '20241124T162727Z_DLLLL_b_0005',
    # '20241124T164605Z_DLLLE_b_0030',
    # '20241124T164731Z_DLLLE_b_0030',
    # '20241124T165036Z_DLLLD_b_0030',
]


@pytest.mark.preprocessed_data
@pytest.mark.parametrize("arm", ["BLUE"])
@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_load_recordings(arm, path_to_legacy_wavecal, path_to_legacy_flats, preprocessed_files_path, etalon_filename):

    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")
    old_flat_file = path_to_legacy_flats / "20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf"

    # Explicit calibration paths (no caldb configured for legacy regression tests)
    calib_root = preprocessed_files_path / 'calibrations'
    flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
    dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

    # read files and instantiate the primitive class
    # this should pick a single file
    raw_file = preprocessed_files_path / (etalon_filename + ".fits")
    # raw_files = sorted([str(f) for f in preprocessed_files_path.glob('20241124T162336Z_DEEEE_*.fits')])
    # selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'WAVECAL', 'BLUE'])

    # Primitives
    adinput = [astrodata.open(raw_file)]

    p = MaroonXSpectrum(adinput)
    p.prepare()
    p.checkArm()
    p.addDQ()  # just placeholder until MX is in caldb
    p.overscanCorrect()
    p.correctImageOrientation()
    p.addVAR(read_noise=True, poisson_noise=True)

    p.extractStripes(flat=flat_path, dark=dark_path)
    adout = p.boxExtraction() # extracts spectra from stripes
    ad = adout[0]

    fibers = None # (1, 2, 3, 4, 5)
    orders = None # (98, 99, 100, 101, 102)
    guess_file = None
    use_sigma_lr = False

    new_generator = load_recordings(ad, guess_file, fibers, orders)
    old_generator = legacy_adapter.load_recordings_legacy(
        str(old_file), str(old_flat_file), guess_file, fibers, orders, use_sigma_lr
        )

    new_dict = dict()
    for fiber, order, data, guess in new_generator:
        new_dict[(int(fiber), int(order))] = (data, guess)

    old_dict = dict()
    for fiber, order, data, guess in old_generator:
        old_dict[(int(fiber), int(order))] = (data, guess)

    # iterate old keys and compare with new keys values
    for (old_fiber, old_order), (old_data, old_guess) in old_dict.items():

        assert (old_fiber, old_order) in new_dict, f"Old fiber/order ({old_fiber}, {old_order}) not found in new data"
        new_data, new_guess = new_dict[(old_fiber, old_order)]

        # Test data shapes are equal
        assert old_data.shape == new_data.shape
        # Test data values are equal
        np.testing.assert_allclose(old_data, new_data, rtol=0, atol=1e-5)


@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_iterative_fit_ETALON(path_to_legacy_wavecal, etalon_filename):
    # Load old peak data. columns are lowercase
    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    old_mask = (old_peak_data["fiber"]==2) & (old_peak_data["order"]==111)

    old_file_inputs = path_to_legacy_wavecal / (etalon_filename + "_2_111_iterative_fit.npy")
    old_input = np.load(old_file_inputs, allow_pickle=True).item()

    set_logger(log)

    # Calculate expected results
    output = maroonx_fit.iterative_fit(
        input_spectrum = old_input['data'],
        degree_sigma = old_input['degree_sigma'],
        degree_width = old_input['degree_width'],
        iterations = old_input['iterations'],
        guess_spectrum = old_input['initial_parameters'],
        fiber = old_input['fiber'],
        plot_path = old_input['plot_path'],
        use_sigma_lr = old_input['use_sigma_lr'],
        show_plots = old_input['show_plots']
        )
    output.fiber = 2
    output.order = 111
    output.recording_time = 0.0
    results = [output]

    peaks = maroonx_fit.insert_peak_parameters(results)

    assert peaks.shape == old_peak_data[old_mask].shape, "Peak data shape mismatch"

    relevant_columns = ['center', 'amplitude', 'offset', 'width', 'sigma1', 'sigma2']

    # Test column values
    for col in relevant_columns:
        np.testing.assert_allclose(
            peaks[col.upper()].values,
            old_peak_data[old_mask][col.lower()].values,
            rtol=0, atol=1e-2)



def test_iterative_fit_LFC(path_to_legacy_wavecal):
    # Load old peak data. columns are lowercase
    old_file = path_to_legacy_wavecal / "20241124T030436Z_DLLLL_b_0005.hdf"
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    old_mask = (old_peak_data["fiber"]==2) & (old_peak_data["order"]==111)

    old_file_inputs = path_to_legacy_wavecal / "20241124T030436Z_DLLLL_b_0005_2_111_iterative_fit.npy"
    old_input = np.load(old_file_inputs, allow_pickle=True).item()

    set_logger(log)

    # Calculate expected results
    output = maroonx_fit.iterative_fit(
        input_spectrum = old_input['data'],
        degree_sigma = old_input['degree_sigma'],
        degree_width = old_input['degree_width'],
        iterations = old_input['iterations'],
        guess_spectrum = old_input['initial_parameters'],
        fiber = old_input['fiber'],
        plot_path = old_input['plot_path'],
        use_sigma_lr = old_input['use_sigma_lr'],
        show_plots = old_input['show_plots']
        )
    output.fiber = 2
    output.order = 111
    output.recording_time = 0.0
    results = [output]

    peaks = maroonx_fit.insert_peak_parameters(results)

    assert peaks.shape == old_peak_data[old_mask].shape, "Peak data shape mismatch"

    relevant_columns = ['center', 'amplitude', 'offset', 'width', 'sigma1', 'sigma2']

    # Test column values
    for col in relevant_columns:
        np.testing.assert_allclose(
            peaks[col.upper()].values,
            old_peak_data[old_mask][col.lower()].values,
            rtol=0, atol=1e-4)  # value * rtol + atol



@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_getPeaksAndPolynomials(path_to_legacy_wavecal, preprocessed_files_path, etalon_filename):

    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    # Load old peak data. columns are lowercase
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    assert "fiber" in old_peak_data.columns

    if USE_CACHE:
        # Use previously saved data on science_dir
        adout = [astrodata.open(str(preprocessed_files_path / (etalon_filename + "_wavecal.fits")))]
    else:
        # Determine arm from filename (_b_ = BLUE, _r_ = RED)
        arm = 'BLUE' if '_b_' in etalon_filename else 'RED'
        calib_root = preprocessed_files_path / 'calibrations'
        flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
        dark_path = None # etalon is not dark corrected

        # Primitives
        adinput = [astrodata.open(str(preprocessed_files_path / (etalon_filename + ".fits")))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        
        p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
        #p.extractStripes(flat=flat_path, dark=dark_path)
        
        p.boxExtraction() # extracts spectra from stripes

        adout = p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))

    # Load new peak data. columns are uppercase
    new_peak_data = adout[0][0].PEAKS.to_pandas()

    # Make sure they are correctly sorted by fiber, order, and center
    new_peak_data.sort_values(by=['FIBER', 'ORDER', 'CENTER'], inplace=True)
    old_peak_data.sort_values(by=['fiber', 'order', 'center'], inplace=True)

    # Test shapes are equal
    assert old_peak_data.shape == new_peak_data.shape

    # Test the the column FIBER has the same value counts
    assert new_peak_data["FIBER"].value_counts().equals(old_peak_data["fiber"].value_counts())

    col_tolerances = {
        'amplitude':  {'atol': 1e-6, 'rtol': 0},
        'center':     {'atol': 1e-6, 'rtol': 0},
        'fiber':      {'atol': 0,    'rtol': 0},
        'order':      {'atol': 0,    'rtol': 0},
        'offset':     {'atol': 1e-6, 'rtol': 0},
        'sigma1':     {'atol': 1e-6, 'rtol': 0},
        'sigma2':     {'atol': 1e-6, 'rtol': 0},
        'width':      {'atol': 1e-6, 'rtol': 0},
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
            print(f"Column {col} mismatch: {err}")
            raise err


@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_fitAndApplyEtalonWls(path_to_legacy_wavecal, preprocessed_files_path, etalon_filename):

    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    # Load old peak data. columns are lowercase
    old_peak_data = pd.read_hdf(old_file, '/peak_data').reset_index()
    assert "fiber" in old_peak_data.columns

    if USE_CACHE:
        # Use previously saved data on science_dir
        adout = [astrodata.open(str(preprocessed_files_path / (etalon_filename + "_wavecal.fits")))]
    else:
        # Determine arm from filename (_b_ = BLUE, _r_ = RED)
        arm = 'BLUE' if '_b_' in etalon_filename else 'RED'
        calib_root = preprocessed_files_path / 'calibrations'
        flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
        dark_path = str(calib_root / _CALIB_FILES[arm]['dark'])

        # Primitives
        adinput = [astrodata.open(str(preprocessed_files_path / (etalon_filename + ".fits")))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
        #p.extractStripes(flat=flat_path, dark=dark_path)

        p.boxExtraction() # extracts spectra from stripes

        p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))
        p.staticWavelengthSolution()
        adout = p.fitAndApplyEtalonWls()

    # Load new peak data. columns are uppercase
    new_peak_data = adout[0][0].PEAK_DATA.to_pandas()

    # Test relevant columns are present
    # Old table has an extra wavelength column, shapes will be different
    new_cols = {c.lower() for c in new_peak_data.columns}
    old_cols = {c.lower() for c in old_peak_data.columns}
    assert new_cols.issubset(old_cols)

    # Test the the column FIBER has the same value counts
    assert new_peak_data["FIBER"].value_counts().equals(old_peak_data["fiber"].value_counts())

    col_tolerances = {
        # Peak fit parameters (from getPeaksAndPolynomials)
        'amplitude':         {'atol': 1e-6, 'rtol': 0},
        'center':            {'atol': 1e-6, 'rtol': 0},
        'fiber':             {'atol': 0,    'rtol': 0},
        'order':             {'atol': 0,    'rtol': 0},
        'offset':            {'atol': 1e-6, 'rtol': 0},
        'sigma1':            {'atol': 1e-6, 'rtol': 0},
        'sigma2':            {'atol': 1e-6, 'rtol': 0},
        'width':             {'atol': 1e-6, 'rtol': 0},
        # Wavelength derived columns (from fitAndApplyEtalonWls)
        'wavelength_by_thar': {'atol': 1e-8, 'rtol': 0},
        'm':                  {'atol': 0,   'rtol': 0},
        'm_fraction':         {'atol': 1e-6, 'rtol': 0},
        'dispersion_mps':     {'atol': 1e-6, 'rtol': 0},
    }
    shared_cols = new_cols & old_cols
    
    for col, tol in col_tolerances.items():
        try:
            np.testing.assert_allclose(
                new_peak_data[col.upper()].values,
                old_peak_data[col.lower()].values,
                **tol,
                err_msg=f"Column '{col}' mismatch",
            )
        except AssertionError as err:
            print(f"Column {col} mismatch: {err}")
            raise err


@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_dynamicWavelengthSolution(path_to_legacy_wavecal, preprocessed_files_path, etalon_filename):

    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    # Load old peak data. columns are lowercase
    legacy_wls = legacy_adapter.load_dict_from_hdf5(str(old_file), 'wavelengths_dynamic/')


    if USE_CACHE:
        # Use previously saved data on science_dir
        adout = [astrodata.open(str(preprocessed_files_path / (etalon_filename + "_wavecal.fits")))]
    else:
        # Determine arm from filename (_b_ = BLUE, _r_ = RED)
        arm = 'BLUE' if '_b_' in etalon_filename else 'RED'
        calib_root = preprocessed_files_path / 'calibrations'
        flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
        dark_path = None # etalon is not dark corrected            

        # Primitives
        adinput = [astrodata.open(str(preprocessed_files_path / (etalon_filename + ".fits")))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
        #p.extractStripes(flat=flat_path, dark=dark_path)

        p.boxExtraction() # extracts spectra from stripes

        p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))
        p.staticWavelengthSolution()
        adout = p.fitAndApplyEtalonWls()

    fail_counter = 0
    for fiber in range(1, 6):
        new_wls = getattr(adout[0][0], f"WLS_DYNAMIC_FIBER_{fiber}")
        if new_wls.size == 1:
            # non reduced fibers (usually fiber 1) have size 1
            continue

        orders = getattr(adout[0][0], f"REDUCED_ORDERS_FIBER_{fiber}")
        orders = orders.astype(int)

        for idx, order in enumerate(orders):
            legacy_order_wls = legacy_wls[f"fiber_{fiber}"][f"{order}"]
            new_order_wls = new_wls[idx, :]

            try:
                np.testing.assert_allclose(legacy_order_wls, new_order_wls, rtol=1e-8)
                print(f'fiber/order : {fiber}/{order} [OK]')
            except AssertionError:
                fail_counter += 1
                print(f'fiber/order : {fiber}/{order} [FAIL]')
    assert fail_counter == 0

@pytest.mark.slow()
@pytest.mark.preprocessed_data
@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_fiber_drifts_ETALON(path_to_legacy_wavecal, preprocessed_files_path, etalon_filename):

    old_file = path_to_legacy_wavecal / (etalon_filename + ".hdf")

    # Load old drift values
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
        # Use previously saved data on science_dir
        adout = [astrodata.open(str(preprocessed_files_path / (etalon_filename + "_wavecal.fits")))]
    else:
        # Determine arm from filename (_b_ = BLUE, _r_ = RED)
        arm = 'BLUE' if '_b_' in etalon_filename else 'RED'
        calib_root = preprocessed_files_path / 'calibrations'
        flat_path = str(calib_root / _CALIB_FILES[arm]['flat'])
        dark_path = None # etalon is not dark corrected

        # Primitives
        adinput = [astrodata.open(str(preprocessed_files_path / (etalon_filename + ".fits")))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)

        p.extractStripes(flat=flat_path, dark_subtraction_skip_fibers=[1, 2, 3, 4, 5])
        #p.extractStripes(flat=flat_path, dark=dark_path)

        p.boxExtraction() # extracts spectra from stripes

        p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))
        p.staticWavelengthSolution()
        adout = p.fitAndApplyEtalonWls()

    # Test fiber drifts
    drifts = adout[0].fiber_drifts()

    blue_precission = 0.1  # m/s
    red_precission = 0.1  # m/s
    if 'BLUE' in adout[0].tags:
        precission = blue_precission
    elif 'RED' in adout[0].tags:
        precission = red_precission
    else:
        raise ValueError("Unknown arm in adout tags. Expected 'BLUE' or 'RED'.")

    for fiber, drift in enumerate(drifts, start=1):
        legacy_drift = legacy_drifts[f"fiber_{fiber}"]
        # if legacy_drift is None, drift should be None too
        if legacy_drift is None:
            assert drift is None, f"Fiber {fiber} drift mismatch: legacy is None, but new is {drift}"
            continue

        # drifts are float values in m/s, test they are within precission m/s
        assert abs(legacy_drift - drift) < precission, f"Fiber {fiber} drift mismatch: {legacy_drift:.2f} vs {drift:.2f} [m/s]"



