import os

import h5py
import numpy as np
import pandas as pd
import tables

from astropy.io import fits


# =====================================================
# HDF5 helper functions
# =====================================================

def load_recordings_legacy(f_data, f_flat, f_guess, fibers, orders, use_sigma_lr):
    """Iterate over the recordings of the spectra and applies the flat data
    Originally, if a guess file was specified, paramters from the guess file would be loaded
    and passed on, but this usage mode was never fully implemented in the etalon_fit.iterative_fit() function
    Instead passing on the data from the guess file

    :param f_data: HDF5 file with the spectrum data
    :param f_flat: HDF5 file with the flat data
    :param f_guess: HDF5 file to use for guessing parameters, may be None
    :param fibers: If set only the given fibers are loaded
    :param orders: If set only the given orders are loaded
    :param use_sigma_lr: Use different polynomials for sigma on left and right
                         flak. Only required when loading guesses
    :returns: iterator of (fiber, order, data)
    :rtype: iterator
    """
    PATH_SPECTRA = '/box_extraction'
    NODE_PARAMETERS = 'etalon_peak_parameters'
    NODE_PEAKS = 'peaks'
    NODE_POLYNOMIALS = 'polynomials'

    open_h5 = tables.open_file

    if f_guess is not None:
        with open_h5(f_data, 'r') as h5f, open_h5(f_flat, 'r') as h5flat, open_h5(f_guess, 'r') as h5guess:
            for node in h5f.walk_nodes(PATH_SPECTRA, 'Array'):
                fiber = int(node._v_parent._v_name.split('_')[-1])
                order = int(node._v_name)

                if fibers and fiber not in fibers:
                    continue
                if orders and order not in orders:
                    continue

                flat = h5flat.get_node(node._v_pathname)
                data = np.array(node) / np.array(flat)
                #guess = guesses.get((fiber, order), None)
                guess = h5guess.get_node(node._v_pathname)
                guess_data = np.array(guess) / np.array(flat)
                guess_data = guess_data /np.nanmedian(guess_data[500:3500]) * np.nanmedian(data[500:3500])
                yield fiber, order, data, guess_data
    else:
        with open_h5(f_data, 'r') as h5f, open_h5(f_flat, 'r') as h5flat:
            for node in h5f.walk_nodes(PATH_SPECTRA, 'Array'):
                fiber = int(node._v_parent._v_name.split('_')[-1])
                order = int(node._v_name)

                if fibers and fiber not in fibers:
                    continue
                if orders and order not in orders:
                    continue

                flat = h5flat.get_node(node._v_pathname)
                data = np.array(node) / np.array(flat)

                yield fiber, order, data, None

def load_mat(name, filepath):
    with h5py.File(filepath, mode="r") as h5f:
        return h5f[name][()]

def load_dict_from_hdf5(filename, path):
    if not path.endswith("/"):
        path += "/"
    if isinstance(filename, str):
        with h5py.File(filename, 'r', libver='latest') as h5file:
            return recursively_load_dict_contents_from_group(h5file, path)
    else:
        return recursively_load_dict_contents_from_group(filename, path)

def recursively_load_dict_contents_from_group(h5file, path):
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans

def read_header_entries(filename, hdr_keys):
    values = []
    with h5py.File(filename, 'r', libver='latest') as h5f:
        for key in hdr_keys:
            values.append(h5f['header'].attrs.get(key))
    return values