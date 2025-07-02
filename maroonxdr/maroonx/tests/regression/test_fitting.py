

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

import maroonx_instruments  # noqa : important to load adclass tags
from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX
from maroonxdr.maroonx.maroonx_utils import load_recordings

# =========================================================
# FILES TO COMPARE
# =========================================================

OLD_BASE_PATH = Path("/home/martin/Projects/MaroonX/legacy/maroonx_base/data2/")
OLD_FILES_PATH = OLD_BASE_PATH / Path("MaroonX_spectra_reduced/20241124")
OLD_FLAT_FILES_PATH = OLD_BASE_PATH / Path("MaroonX_spectra_reduced/Maroonx_masterframes/202411xx/flats")

NEW_FILES_PATH = Path("/home/martin/Projects/MaroonX/MAROONXDR/science_dir")

SCIENCE_DIR = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir')

# =========================================================
# TESTS
# =========================================================


def test_load_recordings():

    old_file = OLD_FILES_PATH / "20241124T162336Z_DEEEE_b_0030.hdf"
    old_flat_file = OLD_FLAT_FILES_PATH / "20241114T19_masterflat_backgroundsubtracted_FFFFF_b_0007.hdf"

    # read files and instantiate the primitive class
    # this should pick a single file
    raw_files = sorted([str(f) for f in SCIENCE_DIR.glob('20241124T162336Z_DEEEE_*.fits')])
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
        np.testing.assert_allclose(old_data, new_data, rtol=1e-4, atol=1e-4)
        

@pytest.mark.parametrize("arm", ["BLUE"])
def test_getPeaksAndPolynomials(arm):

    old_file = OLD_FILES_PATH / "20241124T162336Z_DEEEE_b_0030.hdf"

    # Load old peak data. columns are lowercase
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    assert "fiber" in old_peak_data.columns

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

    adout = p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))
    
    # Load new peak data. columns are uppercase
    new_peak_data = adout[0][0].PEAKS.to_pandas()

    # Test shapes are equal
    assert old_peak_data.shape == new_peak_data.shape

    # Test the the column FIBER has the same value counts
    assert new_peak_data["FIBER"].value_counts().equals(old_peak_data["fiber"].value_counts())


@pytest.mark.parametrize("arm", ["BLUE"])
def test_fitAndApplyEtalonWls(arm):

    old_file = OLD_FILES_PATH / "20241124T162336Z_DEEEE_b_0030.hdf"

    # Load old peak data. columns are lowercase
    old_peak_data = pd.read_hdf(old_file, '/etalon_peak_parameters/peaks')
    assert "fiber" in old_peak_data.columns

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

    p.getPeaksAndPolynomials(fibers=(2, 3, 4, 5))
    p.staticWavelengthSolution()
    adout = p.fitAndApplyEtalonWls()

    # Load new peak data. columns are uppercase
    new_peak_data = adout[0][0].PEAKS.to_pandas()

    # Test shapes are equal
    assert old_peak_data.shape == new_peak_data.shape

    # Test the the column FIBER has the same value counts
    assert new_peak_data["FIBER"].value_counts().equals(old_peak_data["fiber"].value_counts())   


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