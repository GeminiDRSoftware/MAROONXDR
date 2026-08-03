"""Regression tests for the optimalExtraction primitive.

The sparse STRIPES/F_STRIPES/STRIPES_MASKS inputs of optimalExtraction are
deleted after extraction and cannot be persisted to FITS, so the test rebuilds
them live: it runs the reduce recipe chain up to extractStripes on the staged
2D science frame (master flat and synthetic dark passed explicitly), then runs
optimalExtraction restricted to one science fiber to keep the runtime down.
Box-extracted spectra come for free for all traced fibers and are compared
too.

Usage:
    python -m maroonxdr.maroonx.tests.echelle_extraction.test_optimal_extraction --create-inputs
    python -m maroonxdr.maroonx.tests.echelle_extraction.test_optimal_extraction --create-refs
"""

import os

import numpy as np
import pytest
from gempy.utils import logutils
from recipe_system.testing import ref_ad_factory  # noqa - pytest fixture

import astrodata
import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle

# Science frame, master flat, and synthetic dark, per arm.
datasets = [
    ('20250717T144308Z_SOOOE_b_0300.fits',
     '20250701T171553Z_DDDDF_b_0007_DFFFF_flat.fits',
     '20250717T144308Z_SOOOE_b_0300_synth_dark.fits'),
    ('20250717T144308Z_SOOOE_r_0300.fits',
     '20250701T171553Z_DDDDF_r_0002_DFFFF_flat.fits',
     '20250717T144308Z_SOOOE_r_0300_synth_dark.fits'),
]

# Fiber re-run through optimal extraction (the blessed run used [2, 3, 4];
# the extraction is independent per fiber) and fibers traced by the flat.
OPTIMAL_FIBER = 2


# -- Tests ---------------------------------------------------------------------
@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize('filename, flat, dark', datasets)
def test_optimal_extraction(
    filename, flat, dark, path_to_inputs, ref_ad_factory, change_working_dir
):
    """Extract science spectra and test against the ref."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    flat_path = os.path.join(path_to_inputs, flat)
    dark_path = os.path.join(path_to_inputs, dark)

    with change_working_dir():
        logutils.config(file_name=f'log_{filename.replace(".fits", "")}.txt')
        # reduce recipe chain up to the primitive under test
        p = MAROONXEchelle([ad])
        p.prepare()
        p.checkArm()
        p.addDQ()
        p.overscanCorrect()
        p.correctImageOrientation()
        p.addVAR(read_noise=True, poisson_noise=True)
        p.extractStripes(
            flat=flat_path,
            dark=dark_path,
            dark_subtraction_skip_fibers=[5],
            straylight_removal_fibers=[5],
            report=False,
        )
        p.optimalExtraction(
            optimal_extraction_fibers=[OPTIMAL_FIBER], dark=dark_path
        )
        ad_out = p.streams['main'][0]

    ref_ad = ref_ad_factory(filename.replace('.fits', '_reduced.fits'))

    # Optimally extracted spectrum and variance for the re-run fiber
    for ext in (f'OPTIMAL_REDUCED_FIBER_{OPTIMAL_FIBER}',
                f'OPTIMAL_REDUCED_VAR_{OPTIMAL_FIBER}'):
        np.testing.assert_allclose(
            getattr(ad_out[0], ext), getattr(ref_ad[0], ext),
            rtol=1e-6, atol=1e-8, err_msg=ext,
        )

    # Box-extracted spectra for all traced fibers
    for fiber in [2, 3, 4, 5]:
        np.testing.assert_array_equal(
            getattr(ad_out[0], f'REDUCED_ORDERS_FIBER_{fiber}'),
            getattr(ref_ad[0], f'REDUCED_ORDERS_FIBER_{fiber}'),
            err_msg=f'REDUCED_ORDERS_FIBER_{fiber}',
        )
        for prefix in ('BOX_REDUCED_FIBER', 'BOX_REDUCED_VAR',
                       'BOX_REDUCED_FLAT'):
            ext = f'{prefix}_{fiber}'
            out = getattr(ad_out[0], ext)
            ref = getattr(ref_ad[0], ext)
            # The optimal branch stores dead values as NaN, the box-only
            # branch as 0 (the blessed run extracted fibers 2-4 optimally);
            # check the dead columns agree, then compare values.
            np.testing.assert_array_equal(
                np.isnan(ref) | (ref == 0), np.isnan(out) | (out == 0),
                err_msg=ext,
            )
            np.testing.assert_allclose(
                np.nan_to_num(out), np.nan_to_num(ref),
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
    """Copy the science frames and their calibrations into inputs/."""
    import shutil
    from pathlib import Path

    inputs_dir = Path(_module_data_dir()) / 'inputs'
    inputs_dir.mkdir(parents=True, exist_ok=True)
    preprocessed_dir = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files'

    for names in datasets:
        for filename in names:
            shutil.copy2(preprocessed_dir / filename, inputs_dir / filename)
            print(f'Copied to inputs/:\n    {filename}')


def create_refs_recipe():
    """Copy the blessed reduced science products into refs/."""
    import shutil
    from pathlib import Path

    refs_dir = Path(_module_data_dir()) / 'refs'
    refs_dir.mkdir(parents=True, exist_ok=True)
    preprocessed_dir = Path(os.environ['DRAGONS_TEST']) / 'preprocessed_files'

    for filename, _, _ in datasets:
        ref_name = filename.replace('.fits', '_reduced.fits')
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
