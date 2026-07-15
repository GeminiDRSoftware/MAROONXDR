.. _mxcore:

*******************
Core Legacy Modules
*******************

This section documents the core Python packages used by the MAROON-X pipeline
outside of the DRAGONS primitive layer: the echelle-spectrum data model
(``maroonx_echellespectrum``) and the etalon-fit machinery
(``maroonx_fit``). These two modules do the heavy lifting of modelling and
fitting echelle spectra and, in time, they should be refactored as DRAGONS primitives.

.. note::

   These modules are called from the primitives in
   :ref:`primitives` but can also be used directly for offline analysis.

Echelle Spectrum Package
========================

The ``maroonx_echellespectrum`` package defines the class hierarchy used to
represent extracted echelle spectra and their derived products.
``EchelleSpectrum`` is the base class; ``MaroonXSpectrum``, ``EtalonSpectrum`` and
``FlatSpectrum`` are the concrete subclasses used
throughout the reduction.

EchelleSpectrum (base class)
----------------------------

.. automodule:: maroonxdr.maroonx.maroonx_echellespectrum.echellespectrum
   :members:
   :no-inherited-members:

MaroonXSpectrum
---------------

.. automodule:: maroonxdr.maroonx.maroonx_echellespectrum.maroonxspectrum
   :members:
   :no-inherited-members:

EtalonSpectrum
--------------

.. automodule:: maroonxdr.maroonx.maroonx_echellespectrum.etalonspectrum
   :members:
   :no-inherited-members:

FlatSpectrum
------------

.. automodule:: maroonxdr.maroonx.maroonx_echellespectrum.flatspectrum
   :members:
   :no-inherited-members:

Wavelength Solution
-------------------

.. automodule:: maroonxdr.maroonx.maroonx_echellespectrum.wavelengthsolution
   :members:
   :no-inherited-members:

Spectrum Utilities
------------------

.. automodule:: maroonxdr.maroonx.maroonx_echellespectrum.spectrum_utils
   :members:
   :no-inherited-members:

Etalon Fit Package
==================

The ``maroonx_fit`` package contains the fitting routines used to model the
etalon spectrum, extract per-peak parameters, and describe a fit through the
``MaroonXFit`` object.

Fit Driver
----------

.. automodule:: maroonxdr.maroonx.maroonx_fit.maroonx_fit
   :members:
   :no-inherited-members:

MaroonXFit Object
-----------------

.. automodule:: maroonxdr.maroonx.maroonx_fit.maroon_x_fit_object
   :members:
   :no-inherited-members:

Fit Parameters
--------------

.. automodule:: maroonxdr.maroonx.maroonx_fit.maroonx_fit_parameters
   :members:
   :no-inherited-members:

Fit Spectrum
------------

.. automodule:: maroonxdr.maroonx.maroonx_fit.maroonx_fit_spectrum
   :members:
   :no-inherited-members:
