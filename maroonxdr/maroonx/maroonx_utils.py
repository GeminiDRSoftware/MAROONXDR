'''
MAROONX Utils file.  Contains functions used to load data from reference files.
'''
import os
import logging
import astrodata
import numpy as np
from .lookups import siddb, maskdb, refwavelengthdb

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
        for fiber_number in fibers:
            # Access the data from the input file
            if fiber_number == 1:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_1
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_1
                reduced_err = ad[0].BOX_REDUCED_ERR_1
                reduced_flat = ad[0].BOX_REDUCED_FLAT_1
            if fiber_number == 2:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_2
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_2
                reduced_err = ad[0].BOX_REDUCED_ERR_2
                reduced_flat = ad[0].BOX_REDUCED_FLAT_2
            if fiber_number == 3:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_3
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_3
                reduced_err = ad[0].BOX_REDUCED_ERR_3
                reduced_flat = ad[0].BOX_REDUCED_FLAT_3
            if fiber_number == 4:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_4
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_4
                reduced_err = ad[0].BOX_REDUCED_ERR_4
                reduced_flat = ad[0].BOX_REDUCED_FLAT_4
            if fiber_number == 5:
                reduced_orders = ad[0].REDUCED_ORDERS_FIBER_5
                reduced_fiber = ad[0].BOX_REDUCED_FIBER_5
                reduced_err = ad[0].BOX_REDUCED_ERR_5
                reduced_flat = ad[0].BOX_REDUCED_FLAT_5

            i = 0
            for order in reduced_orders:
            # REDUCED_ORDERS_FIBER_X is a list of order keys
                if orders and order not in orders:
                    # Check if any orders were specified and if the current order is one of them
                    i += 1
                    continue
                # Get each individual 4036 pixel row
                data = reduced_fiber[i] / reduced_flat[i] #Normalize according to flat
                yield fiber_number, order, data, None

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
    refwavelength_dir = os.path.join(os.path.dirname(refwavelengthdb.__file__), 'REFWAVELENGTH')
    db_matches = sorted((k, v) for k, v in refwavelengthdb.refwavelength_dict.items()\
    if arm in k) #Check if there is a reference wavelength file for the given arm
    if db_matches:
        refwavelength = db_matches[-1][1]
    else:
        log.warning(f'No reference wavelength file found for {ad.filename}')
        return None
    return refwavelength if refwavelength.startswith(os.path.sep) else \
        os.path.join(refwavelength_dir, refwavelength)
