.. primitives_index.rst

.. _primitives_index:

********************
Primitives Reference
********************

2D Image Processing
===================

.. toctree::
   :maxdepth: 1

   primitive_addDQ
   primitive_validateData
   primitive_standardizeWCS
   primitive_checkArm
   primitive_checkMaster
   primitive_correctImageOrientation
   primitive_addVAR
   primitive_checkND
   primitive_subtractOverscan
   primitive_stackFramesMXCal
   primitive_stackDarks
   primitive_stackFlats
   primitive_findStripes
   primitive_identifyStripes
   primitive_defineFlatStripes
   primitive_removeStrayLight
   primitive_separateFlatStreams
   primitive_combineFlatStreams
   primitive_splitBundle
   primitive_fitDarkCoefficients

Echelle Extraction
==================

.. toctree::
   :maxdepth: 1

   primitive_createSyntheticDark
   primitive_createSyntheticDarkFromCoeffs
   primitive_extractStripes
   primitive_optimalExtraction
   primitive_measureBlaze
   primitive_boxExtraction

Spectrum Processing
===================

.. toctree::
   :maxdepth: 1

   primitive_staticWavelengthSolution
   primitive_getPeaksAndPolynomials
   primitive_fitAndApplyEtalonWls
   primitive_applyWavelengthSolution
   primitive_combineFibers
   primitive_barycentricCorrection
   primitive_displaySpectra
   primitive_separateArmStreams
   primitive_bundleArmStreams

Calibration Database
====================

.. toctree::
   :maxdepth: 1

   primitive_getProcessedDarkCoeff
   primitive_getProcessedWavecal
   primitive_storeProcessedDarkCoeff
   primitive_storeProcessedWavecal
