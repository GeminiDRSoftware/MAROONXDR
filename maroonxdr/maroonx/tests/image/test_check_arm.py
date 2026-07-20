"""Unit tests for the checkArm primitive."""

import logging

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def data():
    """Small frame data."""
    return np.ones((8, 8), dtype=np.float32)


# -- Tests ---------------------------------------------------------------------
def test_checkArm(caplog, data):
    """The first frame sets the arm. frames of the other arm are tossed."""
    caplog.set_level(logging.DEBUG)

    ads = [make_frame('RED', data), make_frame('RED', data), make_frame('BLUE', data)]
    out = MAROONX([]).checkArm(ads)

    assert len(out) == 2
    assert all('RED' in ad.tags for ad in out)
    assert any(
        'Not all frames taken with the same camera arm' in r.message
        for r in caplog.records
    )


def test_checkArm_all_same_arm(caplog, data):
    """All frames of the same arm are kept, with no warning."""
    caplog.set_level(logging.DEBUG)

    ads = [make_frame('RED', data) for _ in range(3)]
    out = MAROONX([]).checkArm(ads)

    assert len(out) == 3
    assert not any(
        'Not all frames taken with the same camera arm' in r.message
        for r in caplog.records
    )


def test_checkArm_undefined_arm(data):
    """A first frame with no defined arm raises an error."""
    undefined = make_frame('RED', data)
    undefined[0].hdr['ARM'] = 'GREEN'

    with pytest.raises(OSError, match='no defined camera arm'):
        MAROONX([]).checkArm([undefined, make_frame('RED', data)])
