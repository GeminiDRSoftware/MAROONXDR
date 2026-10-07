'''
MAROONX Utils file.  Contains functions used to load data from reference files.
'''
import os
import json
import logging
from datetime import datetime

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt

from lmfit import Parameters, Parameter

import astrodata

from .lookups import siddb, maskdb, wavelengthdb

def load_recordings(ad, guess_file, fibers, orders):
    """
    This function iterates over the recordings of a spectra and applies the flat data.
    It also passes on the data from the guess file
    Parameters
    ----------
    ad: AstroData object
        AstroData object with STRIPES, F_STRIPES extensions.
        These are made by the boxExtraction primitive.
    guess_file: str
        filename of the guess file
    fibers: list of ints
        list of fibers to process.
    orders: list
        list of orders to process
        """

    if guess_file is not None:
        # Open as a fits file
        guess = astrodata.open(guess_file)
        if fibers is None:
            fibers = [1,2,3,4,5]
        for fiber_number in fibers:
            # Access the data from the input file by loading the fits extensions into memory
            if fiber_number == 1:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_1
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_1
                reduced_var = ad[0].BOX_REDUCED_VAR_1
                reduced_flat = ad[0].BOX_REDUCED_FLAT_1
                guess_fiber = guess[0].BOX_REDUCED_FIBER_1
                guess_var = guess[0].BOX_REDUCED_VAR_1
            if fiber_number == 2:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_2
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_2
                reduced_var = ad[0].BOX_REDUCED_VAR_2
                reduced_flat = ad[0].BOX_REDUCED_FLAT_2
                guess_fiber = guess[0].BOX_REDUCED_FIBER_2
                guess_var = guess[0].BOX_REDUCED_VAR_2
            if fiber_number == 3:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_3
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_3
                reduced_var = ad[0].BOX_REDUCED_VAR_3
                reduced_flat = ad[0].BOX_REDUCED_FLAT_3
                guess_fiber = guess[0].BOX_REDUCED_FIBER_3
                guess_var = guess[0].BOX_REDUCED_VAR_3
            if fiber_number == 4:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_4
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_4
                reduced_var = ad[0].BOX_REDUCED_VAR_4
                reduced_flat = ad[0].BOX_REDUCED_FLAT_4
                guess_fiber = guess[0].BOX_REDUCED_FIBER_4
                guess_var = guess[0].BOX_REDUCED_VAR_4
            if fiber_number == 5:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_5
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_5
                reduced_var = ad[0].BOX_REDUCED_VAR_5
                reduced_flat = ad[0].BOX_REDUCED_FLAT_5
                guess_fiber = guess[0].BOX_REDUCED_FIBER_5
                guess_var = guess[0].BOX_REDUCED_VAR_5

                for order in reduced_orders:
                    # REDUCED_ORDERS_FIBER_X is a list of order keys
                    if len(orders) > 0 and order not in orders:
                        # Check if any orders were specified and if the current order is one of them
                        continue
                    # Get each individual 4036 pixel row
                    for fiber_row, flat_row, guess_row \
                    in zip(reduced_fiber, reduced_flat, guess_fiber):
                        # Normalize according to flat. avoid invalid value division warning
                        mask = (flat_row != 0) & np.isfinite(flat_row)
                        data = np.divide(fiber_row, flat_row,
                            out=np.full_like(fiber_row, np.nan),
                            where=mask)
                        # TODO: Check with Andreas if we should compute error too
                        guess_data = np.divide(guess_row, flat_row,
                            out=np.full_like(guess_row, np.nan),
                            where=mask)
                        guess_data = guess_data / np.nanmedian(guess_data[500:3500])*np.nanmedian(data[500:3500])

                        #Function operates as a generator function so we use yield
                        yield fiber_number, order, data, guess_data
    else:
        if fibers is None:
            fibers = [1, 2, 3, 4, 5]
        for fiber in fibers:

            # Access the data from the input file
            reduced_fiber = getattr(ad[0], f"BOX_REDUCED_FIBER_{fiber}")
            reduced_orders = getattr(ad[0], f"REDUCED_ORDERS_FIBER_{fiber}")
            reduced_var = getattr(ad[0], f"BOX_REDUCED_VAR_{fiber}")
            reduced_flat = getattr(ad[0], f"BOX_REDUCED_FLAT_{fiber}")
            
            if reduced_fiber.size == 1:
                logging.warning(f"Fiber {fiber} not found in {ad.filename}")
                continue

            if orders is None:
                orders = reduced_orders
            
            for order in orders:
                i = list(reduced_orders.astype(int)).index(int(order))
                # Normalize according to flat. avoid invalid value division warning
                mask = (reduced_flat[i] != 0) & np.isfinite(reduced_flat[i])
                data = np.divide(reduced_fiber[i], reduced_flat[i],
                    out=np.full_like(reduced_fiber[i], np.nan),
                    where=mask)
                yield fiber, order, data, None

def get_sid_filename(ad):
    """
    Gets stripe ID file for input frame.  SID will not be caldb compliant as it is
    instrument specific.

    Returns
     -------
    str/None: Filename of appropriate sid
    """
    log = logging.getLogger(__name__)
    sid_dir = os.path.join(os.path.dirname(siddb.__file__), 'SID')
    sid = siddb.sid_dict.get(ad.camera()) #Check if there is a Stripe ID file for the arm
    if sid is None:
        log.warning(f'No SID found for {ad.filename}')
        return None
    return sid if sid.startswith(os.path.sep) else \
        os.path.join(sid_dir, sid)

def get_bpm_filename(ad):
    """
    Gets the packaged bad pixel mask for input MX science frame.
    Used by addDQ as the fallback when the calibration database has no
    processed BPM for the frame.

    Returns
    -------
    str/None: Filename of the appropriate bpms
    """
    log = logging.getLogger(__name__)
    bpm_dir = os.path.join(os.path.dirname(maskdb.__file__), 'BPM')
    bpm = maskdb.bpm_dict.get(ad.camera()) #Check if there is a BPM for the arm
    if bpm is None:
        log.warning(f'No BPM found for {ad.filename}')
        return None
    return bpm if bpm.startswith(os.path.sep) else \
        os.path.join(bpm_dir, bpm)

def get_refwavelength_filename(ad):
    """
    Gets reference wavelength file for input MX science frame.
    REF wavelength files are not caldb compliant as they are instrument specific.

    Returns
    -------
    str/None: Filename of the appropriate ref wavelength file
    """
    log = logging.getLogger(__name__)
    wavelength_dir = os.path.join(os.path.dirname(wavelengthdb.__file__), 'WLS')
    #Check if there is a reference wavelength file for the arm
    wavelength = wavelengthdb.refwavelength_dict.get(ad.camera())
    if wavelength is None:
        log.warning(f'No reference wavelength file found for {ad.filename}')
        return None
    return wavelength if wavelength.startswith(os.path.sep) else \
        os.path.join(wavelength_dir, wavelength)

def get_statwavelength_filename(ad):
    """
    Gets static wavelength file for input MX science frame.
    Static wavelength files are not caldb compliant as they are instrument specific.

    Returns
    -------
    str/None: Filename of the appropriate wavelength file
    """
    log = logging.getLogger(__name__)
    wavelength_dir = os.path.join(os.path.dirname(wavelengthdb.__file__), 'WLS')
    #Check if there is a static wavelength file for the arm
    wavelength = wavelengthdb.statwavelength_dict.get(ad.camera())
    if wavelength is None:
        log.warning(f'No static wavelength file found for {ad.filename}')
        return None
    return wavelength if wavelength.startswith(os.path.sep) else \
        os.path.join(wavelength_dir, wavelength)

def load_params_from_fits(file, ext_name='PARAMETERS'):
    """
    Load lmfit parameters from a FITS file.

    The file is first retrieved with get_refwavelength_filename(ad) for the
    specific arm.
    
    Parameters
    ----------
    file : str or Path
        Path to the FITS file
    ext_name : str, optional
        Name of the FITS extension containing parameters, default is 'PARAMETERS'
        
    Returns
    -------
    lmfit.Parameters
        The loaded parameters object
    """
    with fits.open(file) as hdul:
        if ext_name not in [hdu.name for hdu in hdul]:
            raise KeyError(f"Extension {ext_name} not found in {file}")
        
        # Get the parameters HDU
        param_hdu = hdul[ext_name]
        
        # Convert HDU data to Table
        param_table = Table.read(param_hdu)
        
        # Create a new Parameters object
        params = Parameters()
        
        # Add each parameter to the Parameters object
        for i in range(len(param_table)):           
            params.add(
                name=param_table['Name'][i], 
                value=float(param_table['Value'][i]),
                min=float(param_table['Min'][i]),
                max=float(param_table['Max'][i]),
                vary=bool(param_table['Vary'][i]),
                )
        
        return params

def load_refwls_from_fits(file, ext_name=None):
    """
    Load wavelength solution for a specific fiber from a FITS file.

    The file is first retrieved with get_refwavelength_filename(ad) for the
    specific arm.

    Parameters
    ----------
    file : str or Path
        Path to the FITS file
    ext_name : str
        Name of the FITS extension with wavelength solution,
        e.g. 'FIBER_2', 'FIBER_3', etc.
        
    Returns
    -------
    dict
        Dictionary containing wavelength solution attributes
    """
        
    with fits.open(file) as hdul:
        # Check if the fiber extension exists
        if ext_name not in [hdu.name for hdu in hdul]:
            raise KeyError(f"Extension {ext_name} not found in {file}")
        
        # Get the extension
        wls_ext = hdul[ext_name]
              
        # Initialize result dictionary
        res = {
            # scalar values
            'max_x': wls_ext.header['MAXX'],
            'poly_deg_x': wls_ext.header['POLYDEGX'],
            'poly_deg_y': wls_ext.header['POLYDEGY'],
            # array values
            'orders': wls_ext.data['ORDERS'].flatten(),
            'weights': wls_ext.data['WEIGHTS'].flatten(),
            'wavelengths': wls_ext.data['WAVELEN'].flatten(),
            'x_norm': wls_ext.data['XNORM'].flatten(),
        }
        
        return res

def load_statwls_from_fits(file, ext_name=None, orders=None):
    """
    Load wavelength solution for a specific fiber from a FITS file.

    The file is first retrieved with get_refwavelength_filename(ad) for the
    specific arm.

    Parameters
    ----------
    file : str or Path
        Path to the FITS file
    ext_name : str
        Name of the FITS extension with wavelength solution,
        e.g. 'FIBER_2', 'FIBER_3', etc.
    orders : list, optional
        List of orders to load. If None, all orders are loaded.

    Returns
    -------
    dict
        Dictionary containing wavelength solution attributes
    """
        
    with fits.open(file) as hdul:
        # Check if the fiber extension exists
        if ext_name not in [hdu.name for hdu in hdul]:
            raise KeyError(f"Extension {ext_name} not found in {file}")
        
        # Get the extension
        wls_ext = hdul[ext_name]
              
        # Initialize result dictionary
        if orders is None:
            all_orders = wls_ext.data.columns.names
            res = {o: wls_ext.data[o] for o in all_orders}
        else:
            orders = [str(int(o)) for o in orders]
            res = {o: wls_ext.data[o] for o in orders}

        return res

def build_bpm_lookup(config_hdf, arm, outdir):
    """
    Build the packaged bad pixel mask lookup from a legacy configuration file.

    The legacy ``/bad_pixel_map`` uses 1 for good pixels; the lookup follows
    the DRAGONS DQ convention (nonzero is bad), so the map is inverted. The
    product is a DRAGONS multi-extension file: primary header plus one
    ``SCI`` extension holding the uint16 mask, tagged as a processed BPM.

    Parameters
    ----------
    config_hdf : str or Path
        Path to the legacy ``config_b.hdf`` or ``config_r.hdf``
    arm : str
        Arm name, ``'BLUE'`` or ``'RED'``
    outdir : str or Path
        Directory where the lookup file is written. The file name is taken
        from ``maskdb.bpm_dict``.

    Returns
    -------
    str
        Path of the written lookup file
    """
    import h5py

    with h5py.File(config_hdf, 'r') as hdf:
        bad_pixel_map = hdf['/bad_pixel_map'][()]
    mask = (1 - bad_pixel_map).astype(np.uint16)

    phu = fits.Header()
    phu['INSTRUME'] = 'MAROON-X'
    phu['OBSTYPE'] = 'BPM'
    phu['ARM'] = arm
    phu['ORIGNAME'] = os.path.basename(config_hdf)
    phu['DATE'] = datetime.now().strftime('%Y-%m-%d')
    phu['PROCBPM'] = (datetime.now().isoformat(timespec='seconds'),
                      'Processed bad pixel mask')

    ad = astrodata.create(phu)
    ad.append(mask, name='SCI', header=fits.Header({'ARM': arm}))

    filename = os.path.join(outdir, maskdb.bpm_dict[arm])
    ad.write(filename, overwrite=True)
    return filename

def build_sid_lookup(config_hdf, arm, outdir):
    """
    Build the stripe identification lookup from a legacy configuration file.

    Writes the ``/identify_stripes`` ``positions`` attribute as a ``SID``
    table with columns ``identify_fiber``, ``fiber_order`` and
    ``fiber_position``.

    Parameters
    ----------
    config_hdf : str or Path
        Path to the legacy ``config_b.hdf`` or ``config_r.hdf``
    arm : str
        Arm name, ``'BLUE'`` or ``'RED'``
    outdir : str or Path
        Directory where the lookup file is written. The file name is taken
        from ``siddb.sid_dict``.

    Returns
    -------
    str
        Path of the written lookup file
    """
    import h5py

    with h5py.File(config_hdf, 'r') as hdf:
        positions = hdf['/identify_stripes'].attrs['positions']

    phu = fits.Header()
    phu['INSTRUME'] = 'MAROON-X'
    phu['OBSTYPE'] = 'SID'
    phu['ARM'] = arm
    phu['ORIGNAME'] = os.path.basename(config_hdf)
    phu['DATE'] = datetime.now().strftime('%Y-%m-%d')

    ad = astrodata.create(phu)
    ad.SID = Table(positions.astype(np.int16),
                   names=('identify_fiber', 'fiber_order', 'fiber_position'))

    filename = os.path.join(outdir, siddb.sid_dict[arm])
    ad.write(filename, overwrite=True)
    return filename

def build_statwls_lookup(config_hdf, arm, outdir):
    """
    Build the static wavelength solution lookup from a legacy configuration
    file.

    Writes one ``FIBER_N`` table per fiber (1 to 5) with one column per
    order, named by the order number, holding the ``/wavelengths_static``
    wavelengths in nm.

    Parameters
    ----------
    config_hdf : str or Path
        Path to the legacy ``config_b.hdf`` or ``config_r.hdf``
    arm : str
        Arm name, ``'BLUE'`` or ``'RED'``
    outdir : str or Path
        Directory where the lookup file is written. The file name is taken
        from ``wavelengthdb.statwavelength_dict``.

    Returns
    -------
    str
        Path of the written lookup file
    """
    import h5py

    phu = fits.Header()
    phu['INSTRUME'] = 'MAROON-X'
    phu['OBSTYPE'] = 'WLSTAT'
    phu['ARM'] = arm
    phu['ORIGNAME'] = os.path.basename(config_hdf)
    phu['DATE'] = datetime.now().strftime('%Y-%m-%d')

    ad = astrodata.create(phu)
    with h5py.File(config_hdf, 'r') as hdf:
        for fiber in [1, 2, 3, 4, 5]:
            orders = hdf[f'/wavelengths_static/fiber_{fiber}']
            setattr(ad, f'FIBER_{fiber}',
                    Table({order: dataset[()] for order, dataset in orders.items()}))

    filename = os.path.join(outdir, wavelengthdb.statwavelength_dict[arm])
    ad.write(filename, overwrite=True)
    return filename

def build_refwls_lookup(peakmodel_hdf, outdir):
    """
    Build the reference wavelength solution lookups, one per arm, from the
    legacy etalon peak model file.

    Each file holds the ``PARAMETERS`` table (the full lmfit Parameter
    record: ``Name``, ``Value``, ``Min``, ``Max``, ``Stderr``, ``Vary``,
    ``Expr``, ``Brute_Step``), which is the same for both arms, and one
    single-row ``FIBER_N`` table per fiber (2 to 4) with the array columns
    ``XNORM``, ``ORDERS``, ``WEIGHTS`` and ``WAVELEN`` and the scalars
    ``MAXX``, ``POLYDEGX`` and ``POLYDEGY`` in the table header.

    Parameters
    ----------
    peakmodel_hdf : str or Path
        Path to the legacy ``wl_combined_final_etalon_peakmodel_2020.hdf``
    outdir : str or Path
        Directory where the lookup files are written. The file names are
        taken from ``wavelengthdb.refwavelength_dict``.

    Returns
    -------
    list of str
        Paths of the written lookup files, blue then red
    """
    import h5py

    filenames = []
    with h5py.File(peakmodel_hdf, 'r') as hdf:
        # The legacy JSON dump also carries dill-encoded callables saved
        # under Python 3.7 that current lmfit cannot decode, so the
        # parameter states are read directly instead of through
        # Parameters.loads. No parameter has an expression, so the
        # callables are not needed.
        parstates = json.loads(hdf['dispersion/parameter'][()])['params']
        params = Parameters()
        for parstate in parstates:
            par = Parameter(name='')
            par.__setstate__(parstate)
            params.add(par)
        params_table = Table({
            'Name': np.array(list(params.keys()), dtype='U20'),
            'Value': [p.value for p in params.values()],
            'Min': [p.min for p in params.values()],
            'Max': [p.max for p in params.values()],
            'Stderr': [np.nan if p.stderr is None else p.stderr
                       for p in params.values()],
            'Vary': [p.vary for p in params.values()],
            'Expr': np.array([str(p.expr) for p in params.values()], dtype='U20'),
            'Brute_Step': [np.nan if p.brute_step is None else p.brute_step
                           for p in params.values()],
        })

        for arm in ['BLUE', 'RED']:
            phu = fits.Header()
            phu['INSTRUME'] = 'MAROON-X'
            phu['OBSTYPE'] = 'WLREF'
            phu['ARM'] = arm
            phu['ORIGNAME'] = os.path.basename(peakmodel_hdf)
            phu['DATE'] = datetime.now().strftime('%Y-%m-%d')

            ad = astrodata.create(phu)
            ad.PARAMETERS = params_table

            for fiber in [2, 3, 4]:
                wls = hdf[f'wls_{arm.lower()}/fiber_{fiber}']
                # Single-row table: each array is one cell
                table = Table({
                    'XNORM': [wls['x_norm'][()]],
                    'ORDERS': [wls['orders'][()]],
                    'WEIGHTS': [wls['weights'][()]],
                    'WAVELEN': [wls['wavelengths'][()]],
                })
                table.meta['header'] = fits.Header({
                    'MAXX': int(wls['maxx'][()]),
                    'POLYDEGX': int(wls['poly_deg_x'][()]),
                    'POLYDEGY': int(wls['poly_deg_y'][()]),
                })
                setattr(ad, f'FIBER_{fiber}', table)

            filename = os.path.join(outdir, wavelengthdb.refwavelength_dict[arm])
            ad.write(filename, overwrite=True)
            filenames.append(filename)

    return filenames

