"""Unit tests for the addDQ primitive."""

import logging

import astrodata
import numpy as np
import pytest

from maroonxdr.maroonx import maroonx_utils
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame

SHAPE = (4400, 4400)


# -- Helpers -------------------------------------------------------------------
def packaged_mask(ad):
    """The packaged lookup mask for the arm of ``ad``, as a DQ array."""
    return astrodata.open(maroonx_utils.get_bpm_filename(ad))[0].data.astype(np.uint16)


def synthetic_bpm(arm, tmp_path):
    """Write a synthetic BPM for ``arm`` with bad pixels and return the path."""
    mask = np.zeros(SHAPE, dtype=np.uint16)
    mask[10, 20] = 1
    mask[100:110, 200] = 1
    bpm = make_frame(arm, mask)
    path = tmp_path / f'bpm_{arm.lower()}.fits'
    bpm.write(str(path))
    return str(path), mask


# -- Tests ---------------------------------------------------------------------
@pytest.mark.preprocessed_data
@pytest.mark.parametrize('arm', ['RED', 'BLUE'])
def test_addDQ_default_uses_packaged_lookup(caplog, arm):
    """With the default static_bpm, the packaged lookup of the arm is applied."""
    caplog.set_level(logging.DEBUG)

    ad = make_frame(arm, np.zeros(SHAPE, dtype=np.float32))
    filename = ad.filename
    expected = packaged_mask(ad)

    out = MAROONX([ad]).addDQ()[0]

    np.testing.assert_array_equal(out[0].mask, expected)
    assert out[0].mask.any()
    assert 'ADDDQ' in out.phu
    assert out.filename == filename
    assert any('using the packaged lookup' in r.message for r in caplog.records)


@pytest.mark.parametrize('arm', ['RED', 'BLUE'])
def test_addDQ_explicit_static_bpm(tmp_path, arm):
    """An explicit static_bpm file is used as given, not the packaged lookup."""
    path, mask = synthetic_bpm(arm, tmp_path)
    ad = make_frame(arm, np.zeros(SHAPE, dtype=np.float32))

    out = MAROONX([ad]).addDQ(static_bpm=path)[0]

    np.testing.assert_array_equal(out[0].mask, mask)
    assert 'ADDDQ' in out.phu


@pytest.mark.parametrize('arm', ['RED', 'BLUE'])
def test_addDQ_no_static_bpm(arm):
    """With static_bpm=None no static mask is applied and the DQ plane is zero."""
    ad = make_frame(arm, np.zeros(SHAPE, dtype=np.float32))

    out = MAROONX([ad]).addDQ(static_bpm=None)[0]

    assert out[0].mask is not None
    assert not out[0].mask.any()
    assert 'ADDDQ' in out.phu
