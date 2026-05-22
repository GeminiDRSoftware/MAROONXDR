"""Tests for variance extension addition primitives."""

import logging
import os
from copy import deepcopy

import astrodata
import numpy as np
import pytest

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

# -- Test datasets -------------------------------------------------------------
# These bundles are needed for debundling into split-arm files
bundles_needed = [
    'N20250721M6125.fits',
]


# -- Tests ---------------------------------------------------------------------
@pytest.mark.parametrize('filename', ['20250721T170049Z_DDDDE_r_0300.fits'])
def test_var_single(caplog, path_to_inputs, filename):
    """Test addVAR primitive creates variance extension correctly.

    Parameters
    ----------
    caplog : fixture
    filename : str

    Returns
    -------
    None
    """
    caplog.set_level(logging.DEBUG)

    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    p = MAROONX([deepcopy(ad)])
    p.prepare()
    out = p.addVAR(read_noise=True, poisson_noise=True)

    assert len(caplog.records) > 0
    assert out[0].variance is not None
    assert any('read noise variance contribution' in r.message for r in caplog.records)

    # check that variance array is not full of zeros
    assert not np.all(out[0].variance == 0)


# -- Create inputs -------------------------------------------------------------
def create_inputs():
    """
    Create input files for this test module.

    Run with: python -m maroonxdr.maroonx.tests.image.test_var --create-inputs

    Reads raw bundles from $DRAGONS_TEST/raw_files/ (populated by the
    download_raws nox session) and runs splitBundle to produce debundled
    single-arm files.
    """
    raw_dir = os.path.join(os.environ['DRAGONS_TEST'], 'raw_files')
    input_path = os.path.join(
        os.environ['DRAGONS_TEST'],
        'maroonxdr',
        'maroonx',
        'image',
        'test_var',
        'inputs',
    )
    os.makedirs(input_path, exist_ok=True)

    for filename in bundles_needed:
        raw_path = os.path.join(raw_dir, filename)
        if not os.path.isfile(raw_path):
            print(f'  Skipping {filename}: not in {raw_dir}')
            continue

        ad = astrodata.open(raw_path)
        p = MAROONX([ad])
        split_ads = p.splitBundle()
        for ad_arm in split_ads:
            ad_arm.write(os.path.join(input_path, ad_arm.filename), overwrite=True)
            print(f'  Wrote {ad_arm.filename} to {input_path}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs()
    else:
        pytest.main()
