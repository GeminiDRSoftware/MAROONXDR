"""Unit tests for the checkND primitive."""

import logging

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def data():
    """Small uniform frame data."""
    return np.ones((8, 8), dtype=np.float32)


# -- Tests ---------------------------------------------------------------------
def test_checkND(caplog, data):
    """All frames sharing the first frame's ND setting are kept, with no warning."""
    caplog.set_level(logging.DEBUG)

    ads = [make_frame('RED', data, nd_position=1.5) for _ in range(3)]
    out = MAROONX([]).checkND(ads)

    assert len(out) == 3
    assert not any('Not all frames have ' in r.message for r in caplog.records)


def test_checkND_mismatched_setting(caplog, data):
    """Frames whose ND setting differs from the first frame's are tossed."""
    caplog.set_level(logging.DEBUG)

    ads = [
        make_frame('RED', data, nd_position=1.5),
        make_frame('RED', data, nd_position=1.5),
        make_frame('RED', data, nd_position=0.5),
    ]
    out = MAROONX([]).checkND(ads)

    assert len(out) == 2
    assert all(ad.filter_orientation()['ND'] == 1.5 for ad in out)
    assert any('Not all frames have ' in r.message for r in caplog.records)


def test_checkND_only_first_unique(caplog, data):
    """A first frame that is the only one with its ND setting raises an error."""
    caplog.set_level(logging.DEBUG)

    ads = [
        make_frame('RED', data, nd_position=1.5),
        make_frame('RED', data, nd_position=0.5),
        make_frame('RED', data, nd_position=0.5),
    ]

    with pytest.raises(OSError):
        MAROONX([]).checkND(ads)

    assert any(
        'Only first frame found, of given, with its simcal ND filter setting'
        in r.message
        for r in caplog.records
    )
