"""Unit tests for the splitBundle primitive."""

import astrodata
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

# -- Fixtures ------------------------------------------------------------------
SCIENCE_FIBER_SETUP = ['Sky', 'Target', 'Target', 'Target', 'Etalon']
DARK_FIBER_SETUP = ['Dark', 'Dark', 'Dark', 'Dark', 'Etalon']


def _make_bundle(fiber_setup, archname):
    """Minimal two-extension MaroonX bundle AstroData object.

    A BLUE and a RED extension, each carrying its own ORIGNAME, is what the
    BUNDLE tag resolves from.
    """
    phu = fits.PrimaryHDU()
    phu.header.set('INSTRUME', 'MAROON-X')
    phu.header.set('DATALAB', 'test')
    phu.header.set('EXPTIME', 300.0)
    phu.header.set('ORIGNAME', archname)
    for number, fiber in enumerate(fiber_setup, start=1):
        phu.header.set(f'FIBER{number}', fiber)

    extensions = []
    for arm in ('BLUE', 'RED'):
        ext = fits.ImageHDU(data=np.ones((32, 32), dtype=np.float32), name='SCI')
        ext.header.set('ARM', arm)
        ext.header.set('EXPTIME', 300.0)
        ext.header.set('ORIGNAME', f'00000000T000000Z_DDDDE_{arm[0].lower()}_0300.fits')
        extensions.append(ext)

    ad = astrodata.create(phu, extensions)
    ad.filename = archname
    return ad


@pytest.fixture()
def ad_bundle():
    """Two-extension DDDDE bundle."""
    return _make_bundle(DARK_FIBER_SETUP, 'N00000000M0000.fits')


@pytest.fixture()
def ad_bundle_exposuremeter():
    """Two-extension SOOOE bundle carrying an EXPOSUREMETER table.

    The RAW + SCI tags reached by the science fiber setup are what trigger the
    TZERO2 / TZERO3 keyword rename in splitBundle.
    """
    ad = _make_bundle(SCIENCE_FIBER_SETUP, 'N00000000M0001.fits')
    exposuremeter = Table({
        'Timestamp': [1.0, 2.0],
        'Flux PC Channel': [10, 20],
        'Flux FRD Channel': [30, 40],
    })
    exposuremeter.meta['header'] = fits.Header({'TZERO2': 7.7, 'TZERO3': 1.0})
    ad.EXPOSUREMETER = exposuremeter
    return ad


# -- Tests ---------------------------------------------------------------------
def test_splitBundle(ad_bundle):
    """Test that splitBundle() correctly splits a MAROON-X bundle into separate arms."""
    ad = ad_bundle

    # Verify this is a bundle with 2 extensions
    assert len(ad) == 2, 'Input should be a bundle with 2 extensions'
    assert 'BUNDLE' in ad.tags, 'Input should have BUNDLE tag'

    # Store original filename
    original_filename = ad.filename

    # Run splitBundle
    p = MAROONX([ad])
    split_ads = p.splitBundle()

    # Should produce 2 outputs (one for each arm)
    assert len(split_ads) == 2, 'splitBundle should produce 2 outputs'

    arm_1, arm_2 = split_ads

    # Each output should have only 1 extension
    assert len(arm_1) == 1, 'should have 1 extension'
    assert len(arm_2) == 1, 'should have 1 extension'

    # Collect ARM values
    assert arm_1.arm()[0] == 'BLUE', 'Extension 1 should be BLUE arm'
    assert arm_2.arm()[0] == 'RED', 'Extension 2 should be RED arm'

    # Verify ARCHNAME references the original bundle
    for arm_ad in split_ads:
        assert (
            arm_ad.phu.get('ARCHNAME') == original_filename
        ), 'ARCHNAME should reference original bundle filename'

    # Verify ORIGNAME is set from extension header
    for arm_ad in split_ads:
        origname = arm_ad.phu.get('ORIGNAME')
        assert origname is not None, 'ORIGNAME should be set'
        assert arm_ad.filename == origname, 'Filename should match ORIGNAME'


def test_splitBundle_exposuremeter(ad_bundle_exposuremeter):
    """Test that splitBundle preserves EXPOSUREMETER table metadata."""
    ad_bundle = ad_bundle_exposuremeter

    tzero2 = ad_bundle.EXPOSUREMETER.meta['header']['TZERO2']
    tzero3 = ad_bundle.EXPOSUREMETER.meta['header']['TZERO3']

    p = MAROONX([ad_bundle])
    out = p.splitBundle()

    for ad in out:
        # Verify the zp values are stored
        assert ad.EXPOSUREMETER.meta['header']['ZP_PC'] == tzero2
        assert ad.EXPOSUREMETER.meta['header']['ZP_FRD'] == tzero3

        # The keywords are renamed, so the original ones are gone
        assert 'TZERO2' not in ad.EXPOSUREMETER.meta['header']
        assert 'TZERO3' not in ad.EXPOSUREMETER.meta['header']
