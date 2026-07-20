"""Unit tests for the separateFlatStreams primitive."""

import logging

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame

FDDDF_FIBER_SETUP = ['Flat lamp', 'Dark', 'Dark', 'Dark', 'Flat lamp']
DDDDF_FIBER_SETUP = ['Dark', 'Dark', 'Dark', 'Dark', 'Flat lamp']
DFFFD_FIBER_SETUP = ['Dark', 'Flat lamp', 'Flat lamp', 'Flat lamp', 'Dark']
SCIENCE_FIBER_SETUP = ['Sky', 'Target', 'Target', 'Target', 'Etalon']


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture()
def data():
    """Small uniform frame data."""
    return np.ones((8, 8), dtype=np.float32)


def make_flat(data, fiber_setup):
    """A flat frame with the given fiber setup."""
    return make_frame('RED', data, fiber_setup=fiber_setup)


# -- Tests ---------------------------------------------------------------------
def test_separateFlatStreams(caplog, data):
    """FDDDF flats go to the main stream, DFFFD flats to the DFFFD_flats stream."""
    caplog.set_level(logging.DEBUG)

    flats = [make_flat(data, FDDDF_FIBER_SETUP) for _ in range(3)]
    flats += [make_flat(data, DFFFD_FIBER_SETUP) for _ in range(2)]

    main = MAROONX(flats).separateFlatStreams()

    assert len(main) == 3
    assert all(ad.fiber_setup() == FDDDF_FIBER_SETUP for ad in main)
    assert not any('Not registered as Flat' in r.message for r in caplog.records)


def test_separateFlatStreams_DFFFD_stream(data):
    """The DFFFD flats are stored in the DFFFD_flats stream."""
    flats = [make_flat(data, FDDDF_FIBER_SETUP) for _ in range(3)]
    flats += [make_flat(data, DFFFD_FIBER_SETUP) for _ in range(2)]

    p = MAROONX(flats)
    p.separateFlatStreams()

    assert len(p.streams['DFFFD_flats']) == 2
    assert all(ad.fiber_setup() == DFFFD_FIBER_SETUP for ad in p.streams['DFFFD_flats'])


def test_separateFlatStreams_DDDDF_to_main(data):
    """DDDDF flats also go to the main stream."""
    flats = [make_flat(data, DDDDF_FIBER_SETUP), make_flat(data, DFFFD_FIBER_SETUP)]

    p = MAROONX(flats)
    main = p.separateFlatStreams()

    assert len(main) == 1
    assert main[0].fiber_setup() == DDDDF_FIBER_SETUP
    assert len(p.streams['DFFFD_flats']) == 1
    assert p.streams['DFFFD_flats'][0].fiber_setup() == DFFFD_FIBER_SETUP


def test_separateFlatStreams_mislabeled(caplog, data):
    """A frame that is not a flat is flagged as mislabeled."""
    caplog.set_level(logging.DEBUG)

    flats = [make_flat(data, FDDDF_FIBER_SETUP), make_flat(data, SCIENCE_FIBER_SETUP)]
    MAROONX(flats).separateFlatStreams()

    assert any('Not registered as Flat' in r.message for r in caplog.records)
