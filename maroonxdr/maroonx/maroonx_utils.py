'''
MAROONX Utils file.  Contains functions used to load data from reference files.
'''
import os
import logging

import numpy as np
from astropy.io import fits
from astropy.table import Table

from lmfit import Parameters

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
                reduced_err = ad[0].BOX_REDUCED_ERR_1
                reduced_flat = ad[0].BOX_REDUCED_FLAT_1
                guess_fiber = guess[0].BOX_REDUCED_FIBER_1
                guess_err = guess[0].BOX_REDUCED_ERR_1
            if fiber_number == 2:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_2
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_2
                reduced_err = ad[0].BOX_REDUCED_ERR_2
                reduced_flat = ad[0].BOX_REDUCED_FLAT_2
                guess_fiber = guess[0].BOX_REDUCED_FIBER_2
                guess_err = guess[0].BOX_REDUCED_ERR_2
            if fiber_number == 3:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_3
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_3
                reduced_err = ad[0].BOX_REDUCED_ERR_3
                reduced_flat = ad[0].BOX_REDUCED_FLAT_3
                guess_fiber = guess[0].BOX_REDUCED_FIBER_3
                guess_err = guess[0].BOX_REDUCED_ERR_3
            if fiber_number == 4:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_4
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_4
                reduced_err = ad[0].BOX_REDUCED_ERR_4
                reduced_flat = ad[0].BOX_REDUCED_FLAT_4
                guess_fiber = guess[0].BOX_REDUCED_FIBER_4
                guess_err = guess[0].BOX_REDUCED_ERR_4
            if fiber_number == 5:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_5
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_5
                reduced_err = ad[0].BOX_REDUCED_ERR_5
                reduced_flat = ad[0].BOX_REDUCED_FLAT_5
                guess_fiber = guess[0].BOX_REDUCED_FIBER_5
                guess_err = guess[0].BOX_REDUCED_ERR_5

                for order in reduced_orders:
                    # REDUCED_ORDERS_FIBER_X is a list of order keys
                    if len(orders) > 0 and order not in orders:
                        # Check if any orders were specified and if the current order is one of them
                        continue
                    # Get each individual 4036 pixel row
                    for fiber_row, flat_row, guess_row \
                    in zip(reduced_fiber, reduced_flat, guess_fiber):
                        data = fiber_row / flat_row #Normalize according to flat
                        # TODO: Check with Andreas if we should compute error too
                        guess_data = guess_row / flat_row # Normalize according to flat
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
            reduced_err = getattr(ad[0], f"BOX_REDUCED_ERR_{fiber}")
            reduced_flat = getattr(ad[0], f"BOX_REDUCED_FLAT_{fiber}")
            
            if reduced_fiber.size == 1:
                logging.warning(f"Fiber {fiber} not found in {ad.filename}")
                continue

            if orders is None:
                orders = reduced_orders
            
            for order in orders:
                i = list(reduced_orders.astype(int)).index(int(order))
                data = reduced_fiber[i] / reduced_flat[i] # Normalize according to flat
                yield fiber, order, data, None

def get_sid_filename(ad):
    """
    Gets stripe ID file for input frame.  SID will not be caldb compliant as it is
    instrument specific.

    Returns
     -------
    str/None: Filename of appropriate sid
    """
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    arm = ('b' if 'BLUE' in ad.tags else 'r') #Get appropriate arm
    sid_dir = os.path.join(os.path.dirname(siddb.__file__), 'SID')
    db_matches = sorted((k, v) for k, v in siddb.sid_dict.items()
                            if arm in k) #Check if there is a Stripe ID file for the given arm
    if db_matches:
        sid = db_matches[-1][1]
    else:
        log.warning(f'No SID found for {ad.filename}')
        return None
    return sid if sid.startswith(os.path.sep) else \
        os.path.join(sid_dir, sid)

def get_bpm_filename(ad):
    """
    Gets bad pixel mask for input MX science frame.
    this function can be removed when MX is bpm caldb compliant

    Returns
    -------
    str/None: Filename of the appropriate bpms
    """
    arm = ('b' if 'BLUE' in ad.tags else 'r') #Get appropriate arm
    bpm_dir = os.path.join(os.path.dirname(maskdb.__file__), 'BPM')
    bpm = 'BPM_'+arm+'_0000.fits' #Append appropriate arm to the bpm name
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
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    arm = ('b' if 'BLUE' in ad.tags else 'r') #Get appropriate arm
    wavelength_dir = os.path.join(os.path.dirname(wavelengthdb.__file__), 'WLS')
    db_matches = sorted((k, v) for k, v in wavelengthdb.refwavelength_dict.items()\
    if arm in k) #Check if there is a reference wavelength file for the given arm
    if db_matches:
        wavelength = db_matches[-1][1]
    else:
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
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    arm = ('b' if 'BLUE' in ad.tags else 'r') #Get appropriate arm
    wavelength_dir = os.path.join(os.path.dirname(wavelengthdb.__file__), 'WLS')
    db_matches = sorted((k, v) for k, v in wavelengthdb.statwavelength_dict.items()\
    if arm in k) #Check if there is a reference wavelength file for the given arm
    if db_matches:
        wavelength = db_matches[-1][1]
    else:
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
            'maxx': wls_ext.header['MAXX'],
            'poly_deg_x': wls_ext.header['POLYDEGX'],
            'poly_deg_y': wls_ext.header['POLYDEGY'],
            # array values
            'orders': wls_ext.data['ORDERS'].flatten(),
            'weights': wls_ext.data['WEIGHTS'].flatten(),
            'wavelengths': wls_ext.data['WAVELEN'].flatten(),
            'x_norm': wls_ext.data['X_NORM'].flatten(),
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