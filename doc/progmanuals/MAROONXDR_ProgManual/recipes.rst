.. _recipes:

*******
Recipes
*******

This section documents all the reduction recipes available in the MAROONX Data Reduction pipeline.
Recipes define the complete processing workflows by combining sequences of primitives to achieve
specific data reduction goals.

.. note::

   Recipes are the primary interface for data reduction. Each recipe is designed for a specific
   observation type or calibration purpose. The recipes are organized by processing mode and
   data type, following the DRAGONS framework conventions.

Science Quality Recipes
========================

The Science Quality (SQ) recipes provide the highest quality data reduction for scientific analysis.

Bundle Processing Recipes
--------------------------

.. automodule:: maroonxdr.maroonx.recipes.sq.recipes_BUNDLE
   :members:
   :undoc-members:
   :show-inheritance:

Calibration Recipes
-------------------

Dark Frame Processing
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: maroonxdr.maroonx.recipes.sq.recipes_DARK
   :members:
   :undoc-members:
   :show-inheritance:

Flat Field Processing
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: maroonxdr.maroonx.recipes.sq.recipes_FLAT_SPECT
   :members:
   :undoc-members:
   :show-inheritance:

Wavelength Calibration Recipes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: maroonxdr.maroonx.recipes.sq.recipes_STATICWAVECAL
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: maroonxdr.maroonx.recipes.sq.recipes_DYNAMIC_WAVECAL
   :members:
   :undoc-members:
   :show-inheritance:

Science Data Recipes
---------------------

.. automodule:: maroonxdr.maroonx.recipes.sq.recipes_ECHELLE_SPECT
   :members:
   :undoc-members:
   :show-inheritance: