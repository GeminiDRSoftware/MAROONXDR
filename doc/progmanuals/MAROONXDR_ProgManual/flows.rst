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

   * QA flows

Bundle Processing Flow
======================

The ``processBundle`` recipe is the entry point for every MAROON-X raw
file. It splits the two-arm GOA bundle into per-arm raw frames that the
rest of the pipeline consumes.

.. graphviz:: flows/processBundle.dot
   :align: center
   :caption: ``processBundle`` recipe flow (SCIENCE bundle shown; the recipe accepts any BUNDLE-tagged input).

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

Dark Coefficient Flow
~~~~~~~~~~~~~~~~~~~~~

The ``makeDarkCoefficients`` recipe fits a per-pixel log-linear model
across a curated selection of processed darks and stores the resulting
coefficient arrays as a calibration. The input is a **selection of
processed darks for a single arm** spanning multiple exposure times,
typically a set of ``[60, 120, 300, 600, 900, 1200, 1800]``.
At least five frames are required; ``checkMaster`` rejects any input
that is not a processed dark and ``checkArm`` rejects any frame from a
different arm. The diagram shows the BLUE arm.

.. graphviz:: flows/makeDarkCoefficients.dot
   :align: center
   :caption: ``makeDarkCoefficients`` recipe flow (BLUE arm shown).

Synthetic Dark Flow
~~~~~~~~~~~~~~~~~~~

A synthetic dark is a ``DDDDE`` dark evaluated per pixel from the
coefficient arrays at an exposure time that was never observed directly.
Two recipes produce one; they differ in where the exposure times come
from and in whether the coefficients are looked up or passed in.

The ``makeSyntheticDark`` recipe reads the exposure times off the raw
per-arm science frames it is given, and pulls the ``DARK_COEFF``
calibration from the calibration database (or takes it from the
``dark_coeff`` parameter of ``createSyntheticDark``). The inputs are
grouped by exposure time and arm and one dark is returned per group, in
place of the inputs, named after the first frame of the group; setting
``individual=True`` returns one dark per input frame. Each product is a
deep copy of that frame with the data replaced and the ``FIBER1`` to
``FIBER5`` keywords rewritten to ``DDDDE``. The diagram shows the BLUE
arm.

.. graphviz:: flows/makeSyntheticDark.dot
   :align: center
   :caption: ``makeSyntheticDark`` recipe flow (BLUE arm shown).

The ``makeSyntheticDarksFromCoeffs`` recipe takes the dark coefficients
calibration itself as its input and reads the exposure times off the
``exptime`` parameter of ``createSyntheticDarkFromCoeffs``, one product
per value. There is no ``prepare`` in this chain and no calibration
database lookup: the input is already a processed calibration. Each
product is a copy of it with the ``COEFF_Z0``, ``COEFF_Z1`` and
``LOGEXPTIME`` extensions dropped, ``EXPTIME`` reset in the primary and
extension headers, and the exposure-time field of ``filename`` and
``ORIGNAME`` rewritten. Dropping the ``COEFF_*`` extensions is what makes
the product tag as ``DARK`` rather than ``DARK_COEFF``. The diagram shows
the BLUE arm.

.. graphviz:: flows/makeSyntheticDarksFromCoeffs.dot
   :align: center
   :caption: ``makeSyntheticDarksFromCoeffs`` recipe flow (BLUE arm shown).

Both routes stamp ``SYNTHETIC_DARK_CREATED`` in the primary header, which
is what ``_tag_dark`` reads to apply ``DARK_SYNTH``.
``CalibrationMAROONX.dark`` searches processed darks in two tiers,
``DARK_SYNTH`` first and master darks second, so a synthetic dark at the
right exposure time wins over a master dark of the same exposure time
taken closer in time. Files tagged ``DARK_COEFF`` are filtered out of the
second tier and are never served as darks.

Flat Field Processing Flow
--------------------------

The ``makeProcessedFlatDFFFF`` recipe stacks and combines a per-arm mix
of ``DDDDF`` and ``DFFFD`` raw flats into a single processed
``FFFFF`` flat with all five fibers illuminated. After a shared
prologue, ``separateFlatStreams`` splits the input into two streams
that run through an identical primitive chain in parallel; the streams
are merged at ``combineFlatStreams`` before the final stripe
identification and CalDB registration. The diagram shows the BLUE arm.

.. graphviz:: flows/makeProcessedFlatDFFFF.dot
   :align: center
   :caption: ``makeProcessedFlatDFFFF`` recipe flow (BLUE arm shown). Double-line arrows indicate primitives run on both streams in parallel.


Dynamic Wavelength Calibration Flow
-----------------------------------

The ``makeDynamicWavecal`` recipe derives a per-arm wavelength solution
from a simultaneous etalon calibration frame. After the standard
prepare/DQ/overscan/VAR prologue, ``extractStripes`` pulls the matching
processed flat and processed dark from the calibration database,
``boxExtraction`` yields per-fiber 1D spectra, and the wavelength chain
(``getPeaksAndPolynomials``, ``staticWavelengthSolution``,
``fitAndApplyEtalonWls``) consumes the shipped ``WLSTAT`` and
``REFWAVELENGTH`` lookups to compute and apply the dynamic solution.
The diagram shows the RED arm.

.. graphviz:: flows/makeDynamicWavecal.dot
   :align: center
   :caption: ``makeDynamicWavecal`` recipe flow (RED arm shown).

Science Data Flows
===================

Echelle Spectroscopy Flow
-------------------------

The ``reduce`` recipe is the main science reduction. A debundled
per-arm raw science exposure goes through the standard
prepare/DQ/overscan/VAR prologue, then ``extractStripes`` pulls the
matching processed flat and processed dark from the calibration
database and ``optimalExtraction`` produces per-fiber 1D spectra
(a second, independent dark fetch happens here). The wavelength chain
fits the etalon peaks of the simultaneous calibration fiber, attaches
the static solution from the ``WLSTAT`` lookup, and applies the
drift-corrected solution using the processed wavecal and the
``REFWAVELENGTH`` lookup. ``combineFibers`` adds the virtual sixth
fiber from the three science fibers, and ``barycentricCorrection``
computes the BERV values before the final ``_reduced`` product is
written to the working directory. Five primitives emit PDF diagnostic
reports by default. The diagram shows the BLUE arm.

.. graphviz:: flows/reduce.dot
   :align: center
   :caption: ``reduce`` recipe flow (BLUE arm shown).

Barycentric Correction Flow
---------------------------

The ``applyBarycentricCorrection`` recipe applies a target-specific barycentric
velocity correction to an already-reduced MAROON-X spectrum. It is run after
the main science reduction workflow, on the ``_reduced`` product.

.. todo::

   Document the barycentric correction flow.

Bundle Export Flow
------------------

The ``exportReducedBundle`` recipe is the converge mirror of
``processBundle``: it takes a pair of per-arm reduced science spectra
(one BLUE, one RED) and recombines them into a single FITS bundle
matching the original GOA filename. You can also pass a list of multiple
pairs; the recipe will match BLUE-RED pairs and output all the reduced bundles.
Pairing uses the ``ARCHNAME`` header keyword, which records the original
pre-split bundle filename.
The bundle is written to the working directory only; it is not
registered with the calibration database.

.. graphviz:: flows/exportReducedBundle.dot
   :align: center
   :caption: ``exportReducedBundle`` recipe flow. The double-line arrow indicates the BLUE and RED streams flow in parallel through the split/merge pair.
