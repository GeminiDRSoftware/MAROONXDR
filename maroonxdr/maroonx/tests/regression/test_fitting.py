

from copy import deepcopy
import numpy as np
import pandas as pd
import tables
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
from maroonxdr.maroonx.maroonx_utils import load_recordings
from maroonxdr.maroonx.maroonx_fit import maroonx_fit, set_logger, get_logger

# =========================================================
# FILES TO COMPARE
# =========================================================

# OLD_BASE_PATH = Path("/home/martin/Projects/MaroonX/legacy/maroonx_base/data2/")
# OLD_FILES_PATH = OLD_BASE_PATH / Path("MaroonX_spectra_reduced/20241124")
# OLD_FLAT_FILES_PATH = OLD_BASE_PATH / Path("MaroonX_spectra_reduced/Maroonx_masterframes/202411xx/flats")

# SCIENCE_DIR = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir')


# Set logger
logutils.config(file_name="test_fitting.log", mode="debug", stomp=True)
log = logutils.get_logger("test_fitting")
log.setLevel("DEBUG")

# =========================================================
# TESTS
# =========================================================
USE_CACHE = False

ETALONS = [
    '20241124T030040Z_DEEEE_r_0004',
    # '20241124T030227Z_DEEEE_r_0004',
    # '20241124T030436Z_DLLLL_r_0004',
    # '20241124T030623Z_DLLLL_r_0004',
    # '20241124T162149Z_DEEEE_r_0004',
    # '20241124T162336Z_DEEEE_r_0004',
    # '20241124T162540Z_DLLLL_r_0004',
    # '20241124T162727Z_DLLLL_r_0004',
    # '20241124T030040Z_DEEEE_b_0030',
    # '20241124T030227Z_DEEEE_b_0030',
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


def test_load_recordings(legacy_reduced_path, legacy_flats_path):

    old_file = legacy_reduced_path / "20241124T162336Z_DEEEE_b_0030.hdf"
    old_flat_file = legacy_flats_path / "20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf"

    # read files and instantiate the primitive class
    # this should pick a single file
    raw_files = sorted([str(f) for f in Path().glob('20241124T162336Z_DEEEE_*.fits')])
    selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'WAVECAL', 'BLUE'])

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
    ad = adout[0]

    fibers = None # (1, 2, 3, 4, 5)
    orders = None # (98, 99, 100, 101, 102)
    guess_file = None
    use_sigma_lr = False

    new_generator = load_recordings(ad, guess_file, fibers, orders)
    old_generator = load_recordings_legacy(str(old_file), str(old_flat_file), guess_file, fibers, orders, use_sigma_lr)

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
        np.testing.assert_allclose(old_data, new_data, rtol=0, atol=1e-4)
    
@pytest.mark.skip(reason="This test is for debuging purposes only.")
def test_iterative_fit_legacy():
    # Load old peak data. columns are lowercase
    old_file = OLD_FILES_PATH / "20241124T162336Z_DEEEE_b_0030.hdf"
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    old_mask = (old_peak_data["fiber"]==2) & (old_peak_data["order"]==111)

    old_file_inputs = OLD_FILES_PATH / "20241124T162336Z_DEEEE_b_0030_2_111_iterative_fit.npy"
    old_input = np.load(old_file_inputs, allow_pickle=True).item()

    # Get the old iterative fit function
    import sys
    import os

    MOD_PATH = "/home/martin/Projects/MaroonX/legacy/maroonx_base/reduce/"
    sys.path.append(MOD_PATH)
    from etalon_fit import iterative_fit as legacy_iterative_fit

    # Calculate expected results
    output = legacy_iterative_fit(
        spectrum_ = old_input['data'],
        degree_sigma = old_input['degree_sigma'],
        degree_width = old_input['degree_width'],
        iterations = old_input['iterations'],
        initial_parameters = old_input['initial_parameters'],
        fiber = old_input['fiber'],
        plot_path = old_input['plot_path'],
        use_sigma_lr = old_input['use_sigma_lr'],
        show_plots = old_input['show_plots']
        )


def test_iterative_fit_ETALON(legacy_reduced_path):
    # Load old peak data. columns are lowercase
    old_file = legacy_reduced_path / "20241124T162336Z_DEEEE_b_0030.hdf"
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    old_mask = (old_peak_data["fiber"]==2) & (old_peak_data["order"]==111)

    old_file_inputs = legacy_reduced_path / "20241124T162336Z_DEEEE_b_0030_2_111_iterative_fit.npy"
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



def test_iterative_fit_LFC(legacy_reduced_path):
    # Load old peak data. columns are lowercase
    old_file = legacy_reduced_path / "20241124T030436Z_DLLLL_b_0005.hdf"
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    old_mask = (old_peak_data["fiber"]==2) & (old_peak_data["order"]==111)

    old_file_inputs = legacy_reduced_path / "20241124T030436Z_DLLLL_b_0005_2_111_iterative_fit.npy"
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



@pytest.mark.slow
@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_getPeaksAndPolynomials(legacy_reduced_path, etalon_filename):

    old_file = legacy_reduced_path / (etalon_filename + ".hdf")

    # Load old peak data. columns are lowercase
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    assert "fiber" in old_peak_data.columns

    if USE_CACHE:
        # Use previously saved data on science_dir
        adout = [astrodata.open((etalon_filename + "_wavecal.fits"))]
    else:
        # Primitives
        adinput = [astrodata.open((etalon_filename + ".fits"))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        
        p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
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

    cols = ['amplitude', 'center', 'fiber', # 'lq_cost', # 'lq_status', 
        'offset', 'order', 'sigma1', 'sigma2', 'width']
    for c in cols:
        try: 
            np.testing.assert_allclose(new_peak_data[c.upper()].values, old_peak_data[c.lower()].values, atol=1e-1)
        except AssertionError as err:
            print(f"Column {c} mismatch: {err}")
            raise err

@pytest.mark.slow
@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_fitAndApplyEtalonWls(legacy_reduced_path, etalon_filename):

    old_file = legacy_reduced_path / (etalon_filename + ".hdf")

    # Load old peak data. columns are lowercase
    old_peak_data = pd.read_hdf(old_file, '/peak_data').reset_index()
    assert "fiber" in old_peak_data.columns

    if USE_CACHE:
        # Use previously saved data on science_dir
        adout = [astrodata.open((etalon_filename + "_wavecal.fits"))]
    else:
        # Primitives
        adinput = [astrodata.open((etalon_filename + ".fits"))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        
        p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
        p.boxExtraction() # extracts spectra from stripes

        p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5), multithreading=True)
        p.staticWavelengthSolution()
        adout = p.fitAndApplyEtalonWls()

    # Load new peak data. columns are uppercase
    new_peak_data = adout[0][0].NEW_PEAKS.to_pandas()

    # Test relevant columns are present
    # Old table has an extra wavelength column, shapes will be different
    new_cols = {c.lower() for c in new_peak_data.columns}
    old_cols = {c.lower() for c in old_peak_data.columns}
    assert new_cols.issubset(old_cols)

    # Test the the column FIBER has the same value counts
    assert new_peak_data["FIBER"].value_counts().equals(old_peak_data["fiber"].value_counts())

@pytest.mark.slow
@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_dynamicWavelengthSolution(legacy_reduced_path, etalon_filename):

    old_file = legacy_reduced_path / (etalon_filename + ".hdf")

    # Load old peak data. columns are lowercase
    legacy_wls = load_dict_from_hdf5(str(old_file), 'wavelengths_dynamic/')


    if USE_CACHE:
        # Use previously saved data on science_dir
        adout = [astrodata.open((etalon_filename + "_wavecal.fits"))]
    else:
        # Primitives
        adinput = [astrodata.open((etalon_filename + ".fits"))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        
        p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
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
                np.testing.assert_allclose(legacy_order_wls, new_order_wls, rtol=1e-4)
                print(f'fiber/order : {fiber}/{order} [OK]')
            except AssertionError as err:
                fail_counter += 1
                print(f'fiber/order : {fiber}/{order} [FAIL]')
    assert fail_counter == 0

@pytest.mark.slow
@pytest.mark.parametrize("etalon_filename", ETALONS)
def test_fiber_drifts_ETALON(legacy_reduced_path, etalon_filename):

    old_file = legacy_reduced_path / (etalon_filename + ".hdf")

    # Load old drift values
    hdr_keys = [
        'Drift_Fiber1', 'Drift_Fiber2', 
        'Drift_Fiber3', 'Drift_Fiber4', 'Instrument_Drift',
    ]
    old_header_drifts = read_header_entries(str(old_file), hdr_keys)

    legacy_drifts = {}
    for fiber, raw_entry in enumerate(old_header_drifts, start=1):
        drift = float(raw_entry.decode().split()[0]) if raw_entry is not None else None
        legacy_drifts[f"fiber_{fiber}"] = drift

    
    if USE_CACHE:
        # Use previously saved data on science_dir
        adout = [astrodata.open((etalon_filename + "_wavecal.fits"))]
    else:
        # Primitives
        adinput = [astrodata.open((etalon_filename + ".fits"))]

        p = MaroonXSpectrum(adinput)
        p.prepare()
        p.checkArm()
        p.addDQ()  # just placeholder until MX is in caldb
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        
        p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
        p.boxExtraction() # extracts spectra from stripes

        p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))
        p.staticWavelengthSolution()
        adout = p.fitAndApplyEtalonWls()

    # Test fiber drifts
    drifts = adout[0].fiber_drifts()

    blue_precission = 0.5  # m/s
    red_precission = 1.3  # m/s
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

# =====================================================
# HDF5 helper functions
# =====================================================


def load_recordings_legacy(f_data, f_flat, f_guess, fibers, orders, use_sigma_lr):
    """Iterate over the recordings of the spectra and applies the flat data
    Originally, if a guess file was specified, paramters from the guess file would be loaded
    and passed on, but this usage mode was never fully implemented in the etalon_fit.iterative_fit() function
    Instead passing on the data from the guess file

    :param f_data: HDF5 file with the spectrum data
    :param f_flat: HDF5 file with the flat data
    :param f_guess: HDF5 file to use for guessing parameters, may be None
    :param fibers: If set only the given fibers are loaded
    :param orders: If set only the given orders are loaded
    :param use_sigma_lr: Use different polynomials for sigma on left and right
                         flak. Only required when loading guesses
    :returns: iterator of (fiber, order, data)
    :rtype: iterator
    """

    PATH_SPECTRA = '/box_extraction'
    NODE_PARAMETERS = 'etalon_peak_parameters'
    NODE_PEAKS = 'peaks'
    NODE_POLYNOMIALS = 'polynomials'

    open_h5 = tables.open_file
    #guesses = load_fit_parameters(f_guess, use_sigma_lr) if f_guess else {}
    # with open_h5(f_data, 'r') as h5f, open_h5(f_flat, 'r') as h5flat:
    if f_guess is not None:
        #with open_h5(f_data, 'r') as h5f, open_h5(f_flat, 'r') as h5flat:
        with open_h5(f_data, 'r') as h5f, open_h5(f_flat, 'r') as h5flat, open_h5(f_guess, 'r') as h5guess:
            for node in h5f.walk_nodes(PATH_SPECTRA, 'Array'):
                fiber = int(node._v_parent._v_name.split('_')[-1])
                order = int(node._v_name)

                if fibers and fiber not in fibers:
                    continue
                if orders and order not in orders:
                    continue

                flat = h5flat.get_node(node._v_pathname)
                data = np.array(node) / np.array(flat)
                #guess = guesses.get((fiber, order), None)
                guess = h5guess.get_node(node._v_pathname)
                guess_data = np.array(guess) / np.array(flat)
                guess_data = guess_data /np.nanmedian(guess_data[500:3500]) * np.nanmedian(data[500:3500])
                yield fiber, order, data, guess_data
    else:
        with open_h5(f_data, 'r') as h5f, open_h5(f_flat, 'r') as h5flat:
            for node in h5f.walk_nodes(PATH_SPECTRA, 'Array'):
                fiber = int(node._v_parent._v_name.split('_')[-1])
                order = int(node._v_name)

                if fibers and fiber not in fibers:
                    continue
                if orders and order not in orders:
                    continue

                flat = h5flat.get_node(node._v_pathname)
                data = np.array(node) / np.array(flat)
                
                yield fiber, order, data, None
                #yield fiber, order, np.array(node) , np.array(flat), None

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

def read_header_entries(filename, hdr_keys):
    values = []
    with h5py.File(filename, 'r', libver='latest') as h5f:
        for key in hdr_keys:
            values.append(h5f['header'].attrs.get(key))
    return values