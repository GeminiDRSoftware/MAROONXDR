.. _flows:

*****
Flows
*****

This section describes the data reduction flows (processing sequences) for different types of MAROONX observations.
Each flow represents a complete processing pathway from raw data to science-ready products.

.. note::

   Data reduction flows define the logical sequence of processing steps required for different
   observation types. This section will be expanded with detailed flow diagrams and descriptions
   as the MAROONX pipeline development progresses.

.. todo::

   Complete the documentation for all MAROONX data reduction flows, including:

   * Calibration flows (darks, flats, wavelength calibrations)
   * Science observation flows
   * Quality assessment flows
   * Flow diagrams

Bundle Processing Flow
======================

.. todo::

   Document the bundle processing flow for MAROONX data, including the splitBundle
   primitive and subsequent processing steps.

Calibration Flows
=================

Dark Frame Processing Flow
--------------------------

The ``makeProcessedDark`` recipe stacks a sequence of debundled per-arm
raw darks into a single processed dark, which is then stored on disk and
registered with the calibration database. The recipe runs independently
for each arm; the diagram below shows the BLUE arm.

.. graphviz:: flows/makeProcessedDark.dot
   :align: center
   :caption: ``makeProcessedDark`` recipe flow (BLUE arm shown).

Flat Field Processing Flow
--------------------------

.. todo::

   Document the flat field creation and processing flow.

Wavelength Calibration Flow
---------------------------

Static Wavelength Calibration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo::

   Document the static wavelength calibration flow using reference solutions.

Dynamic Wavelength Calibration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo::

   Document the dynamic wavelength calibration flow using simultaneous calibration sources.

Science Data Flows
===================

Echelle Spectroscopy Flow
-------------------------

.. todo::

   Document the complete echelle spectroscopy processing flow, including:

   * 2D image processing
   * Fiber spectra extraction
   * Wavelength calibration
   * 1D spectrum creation
   * Barycentric velocity correction
   * Final science product generation

Barycentric Correction Flow
---------------------------

The ``applyBarycentricCorrection`` recipe applies a target-specific barycentric
velocity correction to an already-reduced MAROON-X spectrum. It is run after
the main science reduction workflow, on the ``_reduced`` product.

.. todo::

   Document the barycentric correction flow.