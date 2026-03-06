"""Tests for image orientation correction primitives."""

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
    'N20241114M3295.fits',
    'N20241114M3300.fits',
]


# -- Tests ---------------------------------------------------------------------
@pytest.mark.parametrize("filename", ["20241114T181815Z_DFFFD_r_0002.fits"])
def test_correctImageOrientation_does_not_change_red_frames(caplog, path_to_inputs, filename):
    """Test that orientation does not change for raw red frames.

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
    adtest = p.correctImageOrientation().pop()

    np.testing.assert_allclose(adtest[0].data, ad[0].data)
    assert len(caplog.records) > 0
    assert any("set as red" in r.message for r in caplog.records)


@pytest.mark.parametrize("filename", ["20241114T181959Z_DFFFD_b_0008.fits"])
def test_correctImageOrientation_flips_blue_frames(caplog, path_to_inputs, filename):
    """Test that blue frames are flipped along both axes.

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
    adtest = p.correctImageOrientation().pop()

    np.testing.assert_allclose(adtest[0].data, np.fliplr(np.flipud(ad[0].data)))
    assert len(caplog.records) > 0
    assert any("set as blue" in r.message for r in caplog.records)


# -- Create inputs -------------------------------------------------------------
def create_inputs():
    """
    Create input files for this test module.

    Run with: python -m maroonxdr.maroonx.tests.image.test_image_orientation_corrector --create-inputs
    """
    from astrodata.testing import download_from_archive

    input_path = os.path.join(
        os.environ["DRAGONS_TEST"], "maroonxdr", "maroonx",
        "image", "test_image_orientation_corrector", "inputs"
    )
    os.makedirs(input_path, exist_ok=True)

    # Download bundles and run splitBundle to produce debundled files
    for filename in bundles_needed:
        print(f"  Downloading {filename}")
        raw_path = download_from_archive(filename)
        ad = astrodata.open(raw_path)
        p = MAROONX([ad])
        split_ads = p.splitBundle()
        for ad_arm in split_ads:
            ad_arm.write(os.path.join(input_path, ad_arm.filename), overwrite=True)
            print(f"  Wrote {ad_arm.filename} to {input_path}")


if __name__ == '__main__':
    import sys

    if "--create-inputs" in sys.argv[1:]:
        create_inputs()
    else:
        pytest.main([__file__, '-v'])
