"""
Unit tests for splitBundle primitive.

This is a suite of tests to be run with pytest.
"""
import astrodata
import pytest

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX


def test_split_bundle(download_maroonx_file, request_manifest_file):
    """
    Test that splitBundle() correctly splits a MAROON-X bundle into separate arms.
    """
    # Get first FLAT file from manifest
    bundle_filename = request_manifest_file("FLAT", 0)

    # Download file from archive if not present
    download_maroonx_file(bundle_filename)

    ad = astrodata.open(bundle_filename)

    # Verify this is a bundle with 2 extensions
    assert len(ad) == 2, "Input should be a bundle with 2 extensions"
    assert "BUNDLE" in ad.tags, "Input should have BUNDLE tag"

    # Store original filename
    original_filename = ad.filename

    # Run splitBundle
    p = MAROONX([ad])
    split_ads = p.splitBundle()

    # Should produce 2 outputs (one for each arm)
    assert len(split_ads) == 2, "splitBundle should produce 2 outputs"

    arm_1, arm_2 = split_ads

    # Each output should have only 1 extension
    assert len(arm_1) == 1, "should have 1 extension"
    assert len(arm_2) == 1, "should have 1 extension"

    # Collect ARM values
    assert "BLUE" == arm_1.arm()[0], "Extension 1 should be BLUE arm"
    assert "RED" == arm_2.arm()[0], "Extension 2 should be RED arm"

    # Verify ARCHNAME references the original bundle
    for arm_ad in split_ads:
        assert arm_ad.phu.get("ARCHNAME") == original_filename, \
            "ARCHNAME should reference original bundle filename"

    # Verify ORIGNAME is set from extension header
    for arm_ad in split_ads:
        origname = arm_ad.phu.get("ORIGNAME")
        assert origname is not None, "ORIGNAME should be set"
        assert arm_ad.filename == origname, "Filename should match ORIGNAME"

