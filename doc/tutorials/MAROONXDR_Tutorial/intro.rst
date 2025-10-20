.. intro.rst

.. _maroonx_setup:

************************************
MAROON-X DRP: Setup and Installation
************************************

Introduction
============

This guide walks you through setting up the MAROON-X data reduction pipeline development environment.
The MAROON-X pipeline (``maroonxdr``) is built on top of the DRAGONS framework and provides recipes
and primitives for reducing MAROON-X echelle spectroscopy data.

The setup process uses **nox** (a Python task automation tool) to create a complete development environment with
all necessary dependencies, including the DRAGONS framework, scientific libraries, and the MAROON-X instrument packages.

Prerequisites
=============

Before beginning, ensure you have:

* **Python 3.12** installed on your system
* **Git** for cloning repositories
* **Poetry** for dependency management (will be installed automatically by nox if missing)

.. note:: The MAROON-X pipeline requires Python 3.12 specifically. Using other Python versions may cause compatibility issues.


Installation Steps
==================

Step 1: Clone the MAROONXDR Repository
---------------------------------------

First, clone the MAROON-X DRAGONS repository:

.. code-block:: bash

    cd /path/to/your/projects
    git clone https://github.com/GeminiDRSoftware/MAROONXDR.git
    cd MAROONXDR

Step 2: Install Nox
--------------------

Nox is the task automation tool that manages the development environment setup. Install it using pip:

.. code-block:: bash

    # Install nox globally or in your user directory
    pip install --user nox

    # Or install globally (may require sudo)
    sudo apt install nox

    # Verify installation
    nox --version

Step 3: List Available Nox Sessions
------------------------------------

View all available nox sessions to understand what tasks are available:

.. code-block:: bash

    # List all nox sessions
    nox -l

You should see output similar to:

.. code-block:: text

    Sessions defined in MAROONXDR/noxfile.py:

    - devenv -> Create a development environment.
    - devconda -> Create a conda development environment.
    - unit_tests -> Run unit tests.
    - regression_tests -> Run regression tests.
    - complete_tests -> Run complete workflow tests.
    - docs -> Build documentation using Sphinx.
    - docstyle -> Check docstring style using pydocstyle.

The ``devenv`` session is what we'll use to create the development environment.

Step 4: Create the Development Environment
-------------------------------------------

Run the ``devenv`` nox session to automatically set up the complete development environment:

.. code-block:: bash

    # Create development environment (this may take a few minutes)
    nox -s devenv

.. important:: This command performs the following steps:

    1. Creates a Python 3.12 virtual environment at ``venv/``
    2. Clones the DRAGONS framework into ``DRAGONS/`` directory (if not already present)
    3. Installs DRAGONS in development mode
    4. Installs GeminiCalMgr and GeminiObsDB from GitHub
    5. Installs all Python dependencies from ``pyproject.toml``
    6. Installs ``maroonxdr`` and ``maroonx_instruments`` in editable mode

Once the setup completes, activate the virtual environment:

.. code-block:: bash

    # Activate the virtual environment
    source venv/bin/activate

    # Or use the custom prompt name
    source venv/bin/activate mx_dev

Your shell prompt should now show ``(mx_dev)`` indicating the environment is active.

.. note:: The ``devenv`` session automatically sets the ``MAROONX_DRAGONS_TEST`` environment variable to the project root. If you need to run regression tests against legacy pipeline data, please set ``MAROONX_LEGACY_TEST`` to point to your legacy data directory:

    .. code-block:: bash

        export MAROONX_LEGACY_TEST=/path/to/legacy/maroonx/dataX

    ``dataX/`` should be the data directory that conatins the raw and reduced data:

    .. code-block:: bash

        data10/
        ├── logs
        ├── MaroonX_spectra
        └── MaroonX_spectra_reduced


Step 5: Verify Installation
----------------------------

Verify that the installation was successful:

.. code-block:: bash

    # Check that maroonxdr is installed
    python -c "import maroonxdr; print('maroonxdr installed.')"

    # Check that maroonx_instruments is installed
    python -c "import maroonx_instruments; print('maroonx_instruments installed.')"



Understanding the Development Environment
==========================================

Directory Structure After Setup
--------------------------------

After running ``nox -s devenv``, your directory structure will look like this:

.. code-block:: text

    MAROONXDR/
    ├── DRAGONS/                      # DRAGONS framework (cloned by nox)
    │   ├── geminidr/                 # Core DRAGONS primitives
    │   ├── gempy/                    # DRAGONS utilities
    │   └── astrodata/                # AstroData framework
    ├── maroonxdr/                    # MAROON-X pipeline implementation
    │   └── maroonx/                  # Primitives, recipes, tests
    ├── maroonx_instruments/          # MAROON-X instrument definitions
    │   └── maroonx/                  # AstroData class, tags, descriptors
    ├── venv/                         # Virtual environment (created by nox)
    ├── science_dir/                  # Test data directory (if present)
    ├── doc/                          # Documentation
    ├── noxfile.py                    # Nox configuration (this guide)
    ├── pyproject.toml                # Package metadata and dependencies
    └── README.md                     # Repository README



Alternative: Conda Environment
===============================

If you prefer using Conda instead of virtualenv, use the ``devconda`` session:

.. code-block:: bash

    # Create conda development environment
    nox -s devconda

    # Activate the conda environment
    conda activate mx_devconda


Next Steps
==========

Once your development environment is set up, you can proceed to:

* :ref:`maroonx_reduction` - Learn the Python API for data reduction
* :ref:`maroonx_reduction_cli` - Learn the command-line interface for data reduction

Both tutorials will guide you through reducing MAROON-X data from raw frames to science-ready spectra.