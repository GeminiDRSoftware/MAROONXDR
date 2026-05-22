# Function to load old config files in HDF5 format and convert them to FITS format
# This file should produce:
# BPM_b.fits
# BPM_r.fits
# SID_b.fits
# SID_r.fits
# WLSTAT_b.fits
# WLSTAT_r.fits

from pathlib import Path
import h5py
from astropy.io import fits
import numpy as np
from datetime import datetime


def get_BPM_hdu(hdf, color):
    """
    Extracts the BPM data from the HDF5 file and returns it as a FITS HDU.
    """
    # Get the data
    bpm_data = hdf.get('/bad_pixel_map')[()]

    # create primary hdu
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['DATE'] = datetime.now().strftime('%Y-%m-%d')
    primary_hdu.header['ORIGNAME'] = Path(hdf.filename).name
    primary_hdu.header['INSTRUME'] = 'MAROONX'
    primary_hdu.header['OBSTYPE'] = 'BPM'
    primary_hdu.header['ARM'] = color.upper()

    bpm_hdu = fits.ImageHDU(data=bpm_data, name='BPM')
    return fits.HDUList([primary_hdu, bpm_hdu])


def get_SID_hdu(hdf, color):
    """
    Extracts the SID data from the HDF5 file and returns it as a FITS HDU.
    """
    sid_data = hdf.get('/identify_stripes').attrs['positions']

    # create primary hdu
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['DATE'] = datetime.now().strftime('%Y-%m-%d')
    primary_hdu.header['ORIGNAME'] = Path(hdf.filename).name
    primary_hdu.header['INSTRUME'] = 'MAROONX'
    primary_hdu.header['OBSTYPE'] = 'SID'
    primary_hdu.header['ARM'] = color.upper()

    # Create columns for the FITS table
    col1 = fits.Column(name='identify_fiber', format='I', array=np.array(sid_data[:, 0]))
    col2 = fits.Column(name='fiber_order', format='I', array=np.array(sid_data[:, 1]))
    col3 = fits.Column(name='fiber_position', format='I', array=np.array(sid_data[:, 2]))
    
    # Combine columns into a ColDefs object
    columns = fits.ColDefs([col1, col2, col3])
    
    # Create the binary table HDU
    table_hdu = fits.BinTableHDU.from_columns(columns, name='SID')
    return fits.HDUList([primary_hdu, table_hdu])


def get_WLSTAT_hdu(hdf, color):
    """
    Extracts the WLSTAT data from the HDF5 file and returns it as a FITS HDU.
    """
    # Create primary HDU with metadata
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['DATE'] = datetime.now().strftime('%Y-%m-%d')
    primary_hdu.header['ORIGNAME'] = Path(hdf.filename).name
    primary_hdu.header['INSTRUME'] = 'MAROONX'
    primary_hdu.header['OBSTYPE'] = 'WLSTAT'
    primary_hdu.header['ARM'] = color.upper()

    hdu_list = fits.HDUList([primary_hdu])

    for fiber in [1, 2, 3, 4, 5]:
        cols = []
        # Get the data for each fiber
        orders = hdf.get(f'/wavelengths_static/fiber_{fiber}')
        for order, dataset in orders.items():
            data = np.array(dataset[()])
            cols.append(fits.Column(name=str(order), format='D', array=data))
        
        # Define columns and create table HDU
        columns = fits.ColDefs(cols)
        table_hdu = fits.BinTableHDU.from_columns(columns, name=f'FIBER_{fiber}')
        hdu_list.append(table_hdu)
    return hdu_list



if __name__ == "__main__":
    
    # Path to the HDF5 file
    # This should be replaced with the actual path to your HDF5 file
    MX_BASE = Path("/home/martin/Documentos/Projects/MaroonX/maroonx_base")
    CONFIG_PATH = MX_BASE / Path("data3/MaroonX_spectra_reduced/Maroonx_configfiles/202411xx/")
    
    config_files = {
        'blue': CONFIG_PATH / Path("config_b.hdf"),
        'red': CONFIG_PATH / Path("config_r.hdf")
    }

    for color, filepath in config_files.items():
        with h5py.File(filepath, 'r') as hdf:
            bpm_hdu = get_BPM_hdu(hdf, color)
            sid_hdu = get_SID_hdu(hdf, color)
            wlstat_hdu = get_WLSTAT_hdu(hdf, color)    

        bpm_filename = Path() / f"BPM_{color[0].lower()}.fits"
        bpm_hdu.writeto(bpm_filename, overwrite=True)
        
        sid_filename = Path() / f"SID_{color[0].lower()}.fits"
        sid_hdu.writeto(sid_filename, overwrite=True)

        wlstat_filename = Path() / f"WLSTAT_{color[0].lower()}.fits"
        wlstat_hdu.writeto(wlstat_filename, overwrite=True)
