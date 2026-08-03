# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- ThAr (static) wavelength solutions

## [0.4.0] - 2026-08-03

### Added

- Complete 2D to 1D science reduction recipe (`reduce`): stripe extraction with
  dark subtraction and straylight removal, box and optimal extraction,
  drift-corrected wavelength calibration, fiber combination, and barycentric
  correction from exposure-meter flux-weighted timestamps
- Dynamic wavelength calibration recipe (`makeDynamicWavecal`): etalon peak
  fitting and 30-knot cubic spline solutions, stored as processed wavecal
  calibrations
- Dark coefficient fitting (`makeDarkCoefficients`) and synthetic dark
  generation (`makeSyntheticDark`) for arbitrary exposure times
- Blaze function measurement (`makeBlaze`) per fiber of a processed flat
- Bundle handling: `processBundle` splits GOA bundles into arm files,
  `exportReducedBundle` re-bundles reduced spectra
- Integration with the DRAGONS calibration database (user and local
  databases), including the MAROON-X caltypes `processed_wavecal` and
  `processed_dark_coeff`
- Interactive QA recipes with a Bokeh spectrum viewer (`reduceQA`,
  `makeProcessedFlatQA`)
- nox-based development environment and task automation (`devenv`,
  `devconda`, test, docs and packaging sessions)
- Documentation: Tutorial nad User Manual published on Read the Docs, Programmer
  manual built locally, with recipe and primitive reference pages generated
  from the live docstrings
- GitHub Actions testing workflow
- Reference lookup files (BPM, SID, WLS) distributed as `lookups_files.zip`
  with each release
- BSD 3-Clause LICENSE, matching DRAGONS

### Changed

- Test suite rewritten to DRAGONS conventions: unit and regression tests with
  pytest-dragons inputs/refs layout under `DRAGONS_TEST`, plus preprocess
  scripts that build the shared test data
- Straylight removal reworked; the legacy-compatible variant is preserved for
  regression comparisons
- Recipe docstrings cleaned up and README updated for release

### Fixed

- Optimal extraction was saving the error instead of the variance
- Variance propagation in stacked frames

### Removed

- In-repo legacy regression tests; legacy comparisons now live outside the
  repository
- `science_dir/` deprecated in favor of `$DRAGONS_TEST/raw_files/`

## [0.3.0] - 2023-09-04

### Added

- Code to fit polynomials and peaks to the extracted 1D spectra
- Added tests for this code

### Changed

- Reorganized the code into different class structures
- Changed the tagsets to make more sense

## [0.2.0] - 2023-07-13

### Added

- Added addVAR() method to primitives_maroon.py
- Added test for addVAR() method
- Begun implementation of dynamic wavelength calibration
- Added test recipes

### Changed
- Fixed the tests in find_stripes.
- Cleaned code in the stackFramesMXCal() method
-  Modified .gitignore file to ignore .fits files
## [0.1.2] - 2023-06-14

### Added

- Added "complete" tests, simple scripts that create darks and flats that can be used by an end user to test their
installation.
- Added some documentation about how an end user can use these tests to test their installation
- Added method descriptions for the primitives_maroonx_echelle.py file in README.md
- Finished adding comments to the primitives_maroonx_echelle.py file
### Changed
- Changed the import paths of all tests to be absolute so that they are more robust

### Known Issues
- Currently the tests in find_stripes and 2 tests in checkND fail.  Looking into this.

## [0.1.1] - 2023-06-13

### Added

- Added .pylintrc file to project with edits to disable some pylint warnings and use camelCase for method names
- Added method descriptions for primitives in primitives_maroonx_echelle.py to README.md
- Added some explanation on how to ensure correct operation of provided unit tests to README.md

### Changed

- Formatted log messages in primitives_maroonx.py to use fstring in accordance with PyLint
- Changed import statements in test files and for primitives_maroonx.py to be more robust (previous implementation was not working)

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

