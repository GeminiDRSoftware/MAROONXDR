"""Unit tests for the checkMaster primitive."""

import logging

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame

SCIENCE_FIBER_SETUP = ['Sky', 'Target', 'Target', 'Target', 'Etalon']


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def data():
    """Small uniform frame data."""
    return np.ones((8, 8), dtype=np.float32)


def make_processed_dark(data):
    """A DDDDE dark carrying the PRDKCOEF keyword that marks it PROCESSED."""
    ad = make_frame('RED', data)
    ad.phu['PRDKCOEF'] = 'test'
    return ad


# -- Tests ---------------------------------------------------------------------
def test_checkMaster(data):
    """Processed master frames are kept."""
    ads = [make_processed_dark(data), make_processed_dark(data)]
    out = MAROONX([]).checkMaster(ads)

    assert len(out) == 2
    assert all('PROCESSED' in ad.tags for ad in out)


def test_checkMaster_drop_unprocessed(caplog, data):
    """A frame that is not a processed master is dropped."""
    caplog.set_level(logging.DEBUG)

    ads = [make_processed_dark(data), make_frame('RED', data)]
    out = MAROONX([]).checkMaster(ads)

    assert len(out) == 1
    assert all('PROCESSED' in ad.tags for ad in out)
    assert any('Tossing non-master' in r.message for r in caplog.records)


def test_checkMaster_unknown_type(data):
    """A first frame that is neither a dark nor a flat raises an error."""
    ad = make_frame('RED', data, fiber_setup=SCIENCE_FIBER_SETUP)

    with pytest.raises(ValueError, match='unknown master type'):
        MAROONX([]).checkMaster([ad])
