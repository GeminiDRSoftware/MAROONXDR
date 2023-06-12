# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- Inline comments for some methods in primitives_maroonx.py
- Inline comments for all methods in primitives_maroonx_echelle.py
- addVAR() method
- Dynamic wavelength calibration
- Method decsiptions for methods in primitives_maroonx_echelle.py
- Integration with DRAGONS caldb

## [0.1.0] - 2023-06-12

### Added

- Initial Changelog
- Initial README.md
- simple_test.py added for basic test to see that package runs

### Changed

- Filter size changed from 20x20 to 21x21 to ensure compatibility with photutils
- Corrected stackDarks method to be compatible with DRAGONS 3.1

### Removed

- Removed readme.txt (contents moved to README.md)
- Removed previously added installation instructions (contents moved to README.md)