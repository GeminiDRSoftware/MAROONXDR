"""Unit tests for the boxExtraction primitive.

The fast tests run on synthetic sparse stripes. The
preprocessed_data test rebuilds the makeDynamicWavecal chain up to
extractStripes on a real etalon frame (the sparse STRIPES inputs cannot be
persisted to FITS) and compares the box-extracted spectra against the blessed
wavecal product written by preprocess/wavecal.py, whose BOX_* extensions come
from this same primitive.

Usage:
    python -m maroonxdr.maroonx.tests.echelle_extraction.test_box_extraction --create-inputs
    python -m maroonxdr.maroonx.tests.echelle_extraction.test_box_extraction --create-refs
"""

import os

import numpy as np
import pytest
from gempy.utils import logutils
from recipe_system.testing import ref_ad_factory  # noqa: F401 - pytest fixture
from scipy import sparse

import astrodata
import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

from . import make_echelle_frame

# Box extraction sums each column, so a full stripe of constant value gives
# value * height per column.
SIGNAL = 100.0
FLAT_SIGNAL = 1000.0
HEIGHT = 6
FIBERS = [2, 3, 4]
ORDERS = ['90', '91', '92']

# Raw etalon frame, master flat, and blessed wavecal product, per arm.
datasets = [
    ('20250717T163124Z_DEEEE_b_0010.fits',
     '20250701T171553Z_DDDDF_b_0007_DFFFF_flat.fits',
     '20250717T163124Z_DEEEE_b_0010_wavecal.fits'),
    ('20250717T163124Z_DEEEE_r_0004.fits',
     '20250701T171553Z_DDDDF_r_0002_DFFFF_flat.fits',
     '20250717T163124Z_DEEEE_r_0004_wavecal.fits'),
]

# Fibers traced by the DFFFF master flat.
TRACED_FIBERS = [2, 3, 4, 5]


def make_band(value, height=HEIGHT, ny=10, nx=8):
    """A sparse stripe: a horizontal band of constant value, ``height`` rows tall."""
    array = np.zeros((ny, nx), dtype=np.float32)
    start = (ny - height) // 2
    array[start:start + height, :] = value
    return sparse.csc_matrix(array)


def make_stripes(fibers, orders, value):
    """Build a STRIPES-style dict: {fiber_key: {order_key: sparse band}}."""
    return {
        f'fiber_{fiber}': {order: make_band(value) for order in orders}
        for fiber in fibers
    }


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def ad_stripes():
    """An echelle frame carrying synthetic science, flat, and mask stripes."""
    ad = make_echelle_frame('RED')
    ad[0].STRIPES = make_stripes(FIBERS, ORDERS, SIGNAL)
    ad[0].F_STRIPES = make_stripes(FIBERS, ORDERS, FLAT_SIGNAL)
    ad[0].STRIPES_MASKS = make_stripes(FIBERS, ORDERS, 1.0)
    return ad


# -- Tests ---------------------------------------------------------------------
def test_boxExtraction(ad_stripes):
    """Box-extracted flux equals the column sum of each stripe."""
    result = MAROONXEchelle([]).boxExtraction([ad_stripes]).pop()

    gain = ad_stripes.gain()[0][0]

    for fiber in FIBERS:
        box = getattr(result[0], f'BOX_REDUCED_FIBER_{fiber}')
        flat = getattr(result[0], f'BOX_REDUCED_FLAT_{fiber}')
        var = getattr(result[0], f'BOX_REDUCED_VAR_{fiber}')
        orders = getattr(result[0], f'REDUCED_ORDERS_FIBER_{fiber}')

        assert box.shape == (len(ORDERS), 8)
        np.testing.assert_allclose(box, SIGNAL * HEIGHT)
        np.testing.assert_allclose(flat, FLAT_SIGNAL * HEIGHT)
        np.testing.assert_allclose(var, SIGNAL * HEIGHT / gain)
        np.testing.assert_array_equal(orders, [float(o) for o in ORDERS])


def test_boxExtraction_empty_fibers(ad_stripes):
    """Fibers without stripes get an empty placeholder extension."""
    result = MAROONXEchelle([]).boxExtraction([ad_stripes]).pop()

    for fiber in [1, 5]:
        box = getattr(result[0], f'BOX_REDUCED_FIBER_{fiber}')
        assert box.shape == (1, 1)


@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize('filename, flat, ref_name', datasets)
def test_boxExtraction_real_etalon(
    filename, flat, ref_name, path_to_inputs, ref_ad_factory, change_working_dir
):
    """Box extraction of a real etalon frame is stable against the reference."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    flat_path = os.path.join(path_to_inputs, flat)

    with change_working_dir():
        logutils.config(file_name=f'log_{filename.replace(".fits", "")}.txt')
        # makeDynamicWavecal chain up to the primitive under test
        p = MAROONXEchelle([ad])
        p.prepare()
        p.checkArm()
        p.addDQ()
        p.subtractOverscan()
        p.trimOverscan()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        p.extractStripes(flat=flat_path)
        p.boxExtraction()
        ad_out = p.streams['main'][0]

    ref_ad = ref_ad_factory(ref_name)

    for fiber in TRACED_FIBERS:
        np.testing.assert_array_equal(
            getattr(ad_out[0], f'REDUCED_ORDERS_FIBER_{fiber}'),
            getattr(ref_ad[0], f'REDUCED_ORDERS_FIBER_{fiber}'),
            err_msg=f'REDUCED_ORDERS_FIBER_{fiber}',
        )
        np.testing.assert_array_equal(
            getattr(ad_out[0], f'BPM_FIBER_{fiber}'),
            getattr(ref_ad[0], f'BPM_FIBER_{fiber}'),
            err_msg=f'BPM_FIBER_{fiber}',
        )
        for prefix in ('BOX_REDUCED_FIBER', 'BOX_REDUCED_VAR',
                       'BOX_REDUCED_FLAT'):
            ext = f'{prefix}_{fiber}'
            actual = getattr(ad_out[0], ext)
            expected = getattr(ref_ad[0], ext)
            np.testing.assert_allclose(
                actual, expected,
                rtol=1e-6, atol=1e-8, err_msg=ext,
            )


# -- Recipes to create the inputs and refs -------------------------------------
def _module_data_dir():
    root = os.environ.get('DRAGONS_TEST')
    if root is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')
    module_name = os.path.splitext(os.path.basename(__file__))[0]
    return os.path.join(
        root, 'maroonxdr', 'maroonx', 'echelle_extraction', module_name
    )


def create_inputs_recipe():
    """Copy the raw etalon frames and master flats into inputs/."""
    import shutil
    from pathlib import Path

    inputs_dir = Path(_module_data_dir()) / 'inputs'
    inputs_dir.mkdir(parents=True, exist_ok=True)
    preprocessed_dir = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files'

    for filename, flat, _ in datasets:
        for name in (filename, flat):
            shutil.copy2(preprocessed_dir / name, inputs_dir / name)
            print(f'Copied to inputs/:\n    {name}')


def create_refs_recipe():
    """Copy the blessed wavecal products into refs/."""
    import shutil
    from pathlib import Path

    refs_dir = Path(_module_data_dir()) / 'refs'
    refs_dir.mkdir(parents=True, exist_ok=True)
    preprocessed_dir = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files'

    for _, _, ref_name in datasets:
        shutil.copy2(preprocessed_dir / ref_name, refs_dir / ref_name)
        print(f'Copied to refs/:\n    {ref_name}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs_recipe()
    elif '--create-refs' in sys.argv[1:]:
        create_refs_recipe()
    else:
        pytest.main([__file__])
