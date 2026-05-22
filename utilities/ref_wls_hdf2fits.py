from datetime import datetime
from pathlib import Path 

from astropy.io import fits
from astropy.table import Table
import h5py
import numpy as np
import pandas as pd

from lmfit import Parameters, Parameter

# Load old parametes
def load_params_from_hdf(file, path='dispersion/parameter'):
    with h5py.File(file, "r") as h5f:
        if path not in h5f.keys():
            raise KeyError(f"Data at {path} not found in {file}")

        else:
            res = Parameters()
            return res.loads(h5f[path][()])
        
# Load old wls
def load_wls_from_hdf(file, path):
    """

    Args:
        file:
        path:

    Returns:

    """
    with h5py.File(file, "r") as h5f:
        if path not in h5f.keys():
            raise KeyError(f"Data at {path} not found in {file}")

        else:
            res = {}
            for attr in ["maxx", "x_norm", "orders", "weights", "wavelengths", "poly_deg_x", "poly_deg_y"]:
                res[attr] = h5f[path+"/"+attr][()]
            return res #WavelengthSolution(**res)
        


def params_to_HDU(params, ext_name='PARAMETERS'):
    """
    Convert lmfit.Parameters object to a FITS binary table HDU.
    
    Parameters
    ----------
    params : lmfit.Parameters
        The parameters object to convert
    ext_name : str, optional
        Name for the FITS extension, default is 'PARAMETERS'
        
    Returns
    -------
    astropy.io.fits.BinTableHDU
        FITS binary table HDU containing the parameter information
    """
    # Create lists to store each column
    names = []
    values = []
    mins = []
    maxs = []
    varies = []
    
    # Loop through each parameter
    for name, param in params.items():
        names.append(name)
        values.append(param.value)
        mins.append(float('-inf') if param.min is None else param.min)
        maxs.append(float('inf') if param.max is None else param.max)
        varies.append(param.vary)

    
    # Create columns for the FITS table
    col1 = fits.Column(name='Name', format='A20', array=np.array(names))
    col2 = fits.Column(name='Value', format='D', array=np.array(values))
    col3 = fits.Column(name='Min', format='D', array=np.array(mins))
    col4 = fits.Column(name='Max', format='D', array=np.array(maxs))
    col5 = fits.Column(name='Vary', format='L', array=np.array(varies))
    
    # Combine columns into a ColDefs object
    columns = fits.ColDefs([col1, col2, col3, col4, col5])
    
    # Create the binary table HDU
    table_hdu = fits.BinTableHDU.from_columns(columns, name=ext_name)
    
    return table_hdu

def wls_to_HDU(wls, ext_name='WLS'):
    """
    Convert wavelength solution attributes to a FITS binary table HDU.
    
    Parameters
    ----------
    wls : dict
        Dictionary containing wavelength solution attributes
    ext_name : str, optional
        Name for the FITS extension, default is 'WLS'
        
    Returns
    -------
    astropy.io.fits.BinTableHDU
        FITS binary table HDU containing the wavelength solution data
    """
    # Separate arrays from scalars
    arrays = {}
    scalars = {}
    
    for attr, value in wls.items():
        if isinstance(value, np.ndarray) and value.size > 1:
            arrays[attr] = value
        else:
            scalars[attr] = value
    
    # Create columns for arrays
    table = Table()
    for name, array in arrays.items():
        col_name = name.upper().replace('_', '')
        if name == 'wavelengths':
            col_name = 'WAVELEN'  # Shorter name for column
            
        # Add as column to table
        table[col_name] = [array]  # Wrap in list to create a single row table
    
    # Convert table to HDU
    table_hdu = fits.BinTableHDU(table, name=ext_name)
    
    # Add scalar values as header keywords
    for name, value in scalars.items():
        keyword = name.upper().replace('_', '')
        table_hdu.header[keyword] = value
    
    return table_hdu

def convert_hdf5_to_hdus(hdf5_file):
    """
    Convert wavelength solution HDF5 file to FITS HDULists.
    
    Args:
        hdf5_file (str): Path to the input HDF5 file
    
    Returns:
        dict: Dictionary with two keys 'BLUE' and 'RED', each containing an HDUList
    """
    # Initialize result dictionary
    result = {}
    
    params = load_params_from_hdf(hdf5_file, path='dispersion/parameter')
    
    # Process each color
    for color_key in ['BLUE', 'RED']:
        # Create primary HDU with metadata
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['DATE'] = datetime.now().strftime('%Y-%m-%d')
        primary_hdu.header['ORIGNAME'] = Path(hdf5_file).name
        primary_hdu.header['COLOR'] = color_key
        
        # Create HDU list with primary HDU
        hdul = fits.HDUList([primary_hdu])
        #result[color_key] = hdul
    
        # Append Parameters HDU
        params_hdu = params_to_HDU(params, ext_name='PARAMETERS')
        hdul.append(params_hdu)

        for fiber in [2, 3, 4]:
            wls = load_wls_from_hdf(hdf5_file, path=f'wls_{color_key.lower()}/fiber_{fiber}')
            wls_hdu = wls_to_HDU(wls, ext_name=f'FIBER_{fiber}')
            hdul.append(wls_hdu)

        # Add the HDU list to the result dictionary
        result[color_key] = hdul

    return result

# NEW functions to save and load parameters and wavelength solutions

def load_params_from_fits(file, ext_name='PARAMETERS'):
    """
    Load lmfit parameters from a FITS file.
    
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

def load_wls_from_fits(file, ext_name=None):
    """
    Load wavelength solution for a specific fiber from a FITS file.
    
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
            'x_norm': wls_ext.data['XNORM'].flatten(),
        }
        
        return res



if __name__ == "__main__":
    
    # Path to the HDF5 file
    # This should be replaced with the actual path to your HDF5 file
    MX_BASE = Path("/home/martin/Documentos/Projects/MaroonX/maroonx_base")
    WL_PATH = MX_BASE / Path("data3/MaroonX_spectra_reduced/Maroonx_wls/202005xx/")
    WL_MODEL = WL_PATH / Path("wl_combined_final_etalon_peakmodel_2020.hdf")

    # Convert the HDF5 file to a dict of HDULists
    hdu_dict = convert_hdf5_to_hdus(WL_MODEL)
    for color, hdul in hdu_dict.items():
        # Save each HDUList to a FITS file
        # REFWAVELENGTH_r.fits, REFWAVELENGTH_b.fits
        output_filename = Path() / f"REFWAVELENGTH_{color[0].lower()}.fits"
        hdul.writeto(output_filename, overwrite=True)

        print(f"Saved {color} wavelength solution to {output_filename}")


def save_params_to_fits(lmfit_parameters, blue_file, red_file, ext_name='PARAMETERS', overwrite=False):
    """
    Save lmfit parameters to both BLUE and RED FITS files.
    
    Parameters
    ----------
    lmfit_parameters : lmfit.Parameters
        The parameters to save
    blue_file : str
        Path to the BLUE FITS file
    red_file : str
        Path to the RED FITS file
    ext_name : str, optional
        Name for the FITS extension, default is 'PARAMETERS'
    overwrite : bool, optional
        Whether to overwrite existing files, default is False
    """
    # Create parameters table HDU
    params_hdu = params_to_fits_table(lmfit_parameters, ext_name)
    
    # Update both BLUE and RED files
    for file_path in [blue_file, red_file]:
        try:
            # Try to open the existing file
            with fits.open(file_path, mode='update') as hdul:
                # Check if the extension already exists and remove it
                if ext_name in [hdu.name for hdu in hdul]:
                    del hdul[ext_name]
                
                # Add the parameters HDU
                hdul.append(params_hdu)
                hdul.flush()
        except FileNotFoundError:
            # If the file doesn't exist, create it with a primary HDU
            primary_hdu = fits.PrimaryHDU()
            primary_hdu.header['CONTENT'] = 'Wavelength Solution'
            color = 'BLUE' if 'blue' in file_path.lower() or '_b.' in file_path.lower() else 'RED'
            primary_hdu.header['COLOR'] = color
            
            # Create HDU list with primary HDU and parameters
            hdul = fits.HDUList([primary_hdu, params_hdu])
            hdul.writeto(file_path, overwrite=overwrite)


def save_wls_to_fits(wls, file, ext_name='WLS', overwrite=False):
    """
    Save wavelength solution to a FITS file.
    
    Parameters
    ----------
    wls : dict
        Dictionary containing wavelength solution attributes
    file : str
        Path to the FITS file
    ext_name : str, optional
        Name for the FITS extension, default is 'WLS'
    overwrite : bool, optional
        Whether to overwrite existing extension, default is False
    """
    # Create wavelength solution table HDU
    wls_hdu = wls_to_fits_extension(wls, ext_name)
    
    try:
        # Try to open the existing file
        with fits.open(file, mode='update') as hdul:
            # Check if the extension already exists and remove it
            if ext_name in [hdu.name for hdu in hdul]:
                del hdul[ext_name]
            
            # Add the wavelength solution HDU
            hdul.append(wls_hdu)
            hdul.flush()
    except FileNotFoundError:
        # If the file doesn't exist, create it with a primary HDU
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['CONTENT'] = 'Wavelength Solution'
        color = 'BLUE' if 'blue' in file.lower() or '_b.' in file.lower() else 'RED'
        primary_hdu.header['COLOR'] = color
        
        # Create HDU list with primary HDU and wavelength solution
        hdul = fits.HDUList([primary_hdu, wls_hdu])
        hdul.writeto(file, overwrite=overwrite)


