import numpy as np
import astrodata
import pandas as pd
from . import maroonx_etalon_fit

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

def insert_peak_parameters(results):
    """
    Extract the parameters for each peak.

    Parameters
    ----------
    results: FitResults NamedTuple
        Results of the fits generated by iterative fit
    
    Returns
    -------
    peaks: A pandas dataframe that can be converted to a fits table
    """
    peaks = pd.DataFrame()
    
    # Making the empty columns
    fiber_list = []
    order_list = []
    offset_list = []
    center_list = []
    amplitude_list = []
    sigma1_list = []
    sigma2_list = []
    width_list = []
    recorded_list = []
    lq_cost_list = []
    lq_status_list = []

    for result in results:
        offset = maroonx_etalon_fit.get_offset(result.parameters, result.meta_parameters)
        fits = maroonx_etalon_fit.eval_polynomials_at_centers(result.parameters, result.meta_parameters)

        for values, f in zip(fits.T, result.peak_fit_results):
            center, amplitude, sigma1, sigma2, width = values
            
            # Append the values to the columns
            fiber_list.append(result.fiber)
            order_list.append(result.order)
            offset_list.append(offset)
            center_list.append(center)
            amplitude_list.append(amplitude)
            sigma1_list.append(sigma1)
            sigma2_list.append(sigma2)
            width_list.append(width)
            recorded_list.append(result.recording_time)
            lq_cost_list.append(f.cost)
            lq_status_list.append(f.status)

    # Insert the columns into the dataframe
    peaks['fiber'] = fiber_list
    peaks['order'] = order_list
    peaks['offset'] = offset_list
    peaks['center'] = center_list
    peaks['amplitude'] = amplitude_list
    peaks['sigma1'] = sigma1_list
    peaks['sigma2'] = sigma2_list
    peaks['width'] = width_list
    peaks['recorded'] = recorded_list
    peaks['lq_cost'] = lq_cost_list
    peaks['lq_status'] = lq_status_list

    return peaks

def insert_polynomial_parameters(results):
    """
    Extract polynomial parameters from results.

    Parameters
    ----------
    results: FitResults NamedTuple
        Results of the fits generated by iterative fit
    
    """
    meta_parameters = results[0].meta_parameters
    
    #Check valid meta_parameters.  Do not check total and number of peaks.
    if not all(result.meta_parameters.offset == meta_parameters.offset for result in results)\
            and all(result.meta_parameters.sigma == meta_parameters.sigma for result in results)\
            and all(result.meta_parameters.width == meta_parameters.width for result in results)\
            and all(result.meta_parameters.use_sigma_lr == meta_parameters.use_sigma_lr for result in results)\
            and all(result.meta_parameters.use_sigma_lr == meta_parameters.use_sigma_lr for result in results):
        raise ValueError('All meta_parameters must be the same')
    
    if meta_parameters.use_sigma_lr == True:
        poly = pd.DataFrame()
    else:
        poly = pd.DataFrame()
    # Make the empty columns
    recorded_list = []
    fiber_list = []
    order_list = []
    fitrange_list = []
    offset_list = []
    lq_cost_list = []
    lq_status_list = []
    if meta_parameters.use_sigma_lr == True:
        sigma_l_coefficients_list = []
        sigma_r_coefficients_list = []
    else:
        sigma_coefficients_list = []
    width_coefficients_list = []

    # Get the polynomial parameters and insert into dataframe
    for result in results:
        values = maroonx_etalon_fit.split_parameters(result.parameters, result.meta_parameters)
        offset, p_sigma_l, p_sigma_r, p_width, _, _ = values

        f = result.polynomial_fit_result

        # Append the values to the columns
        recorded_list.append(result.recording_time)
        fiber_list.append(result.fiber)
        order_list.append(result.order)

        # Convert the fitrange to a string to store in the fits table
        fitrange = [result.x[0], result.x[-1] + 1]
        fitrange_str = ','.join([str(x) for x in fitrange])
        fitrange_list.append(fitrange_str)

        offset_list.append(offset)
        lq_cost_list.append(f.cost)
        lq_status_list.append(f.status)
        if meta_parameters.use_sigma_lr == True:
            # Convert the list of coefficients to a string to store in the fits table
            sigma_l_coefficients_str = ','.join([str(x) for x in list(p_sigma_l)])
            sigma_r_coefficients_str = ','.join([str(x) for x in list(p_sigma_r)])
            sigma_l_coefficients_list.append(sigma_l_coefficients_str)
            sigma_r_coefficients_list.append(sigma_r_coefficients_str)

        else:
            sigma_coefficients_str = ','.join([str(x) for x in list(p_sigma_l)])
            sigma_coefficients_list.append(sigma_coefficients_str)
    
        width_coefficients_str = ','.join([str(x) for x in list(p_width)])
        width_coefficients_list.append(width_coefficients_str)

    # Insert the columns into the dataframe
    poly['recorded'] = recorded_list
    poly['fiber'] = fiber_list
    poly['order'] = order_list
    poly['fitrange'] = fitrange_list
    poly['offset'] = offset_list
    poly['lq_cost'] = lq_cost_list
    poly['lq_status'] = lq_status_list
    if meta_parameters.use_sigma_lr == True:
        poly['sigma_l_coefficients'] = sigma_l_coefficients_list
        poly['sigma_r_coefficients'] = sigma_r_coefficients_list
    else:
        poly['sigma_coefficients'] = sigma_coefficients_list
    poly['width_coefficients'] = width_coefficients_list

    return poly


        