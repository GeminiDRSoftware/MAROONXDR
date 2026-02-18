"""Diagnostic plotting functions for the MaroonX DRAGONS pipeline."""

import numpy as np
import matplotlib.pyplot as plt


def plot_backgroundfit(data_original, data_masked, bkg1, bkg2, result,
                       max_val_masked, bkg1_mesh_nmasked=None, box_size=None):
    """
    Create diagnostic figures for straylight removal.

    Parameters
    ----------
    data_original : ndarray
        Original frame data before straylight removal.
    data_masked : ndarray
        Frame data with orders masked for background fitting.
    bkg1 : Background2D or ndarray
        First iteration background model. If a Background2D object, the
        background array and mesh info are extracted from it. If a plain
        ndarray (legacy path), ``bkg1_mesh_nmasked`` and ``box_size`` must
        be supplied.
    bkg2 : Background2D or ndarray
        Second iteration (correction) background model.
    result : ndarray
        Final straylight-subtracted data.
    max_val_masked : float
        Maximum value in the masked data used for scaling plots.
    bkg1_mesh_nmasked : ndarray, optional
        Number of masked pixels per mesh cell (legacy path only).
    box_size : int, optional
        Background mesh box size in pixels (legacy path only).

    Returns
    -------
    list of matplotlib.figure.Figure
        Seven figures, one per diagnostic stage.
    """
    if hasattr(bkg1, 'background'):
        bkg1_background = bkg1.background
        bkg1_mesh_nmasked = bkg1.mesh_nmasked
        bkg2_background = bkg2.background
        box_size = bkg1.box_size[0]
    else:
        bkg1_background = bkg1
        bkg2_background = bkg2

    figsize = (10, 8)
    figures = []

    # Figure 1: Raw frame, bias corrected
    fig = plt.figure(figsize=figsize)
    plt.title(f'Raw frame, bias corrected (max signal level: {int(max_val_masked):5d})')
    plt.imshow(data_original, origin='lower', vmin=0, vmax=np.nanmax(data_masked))
    plt.colorbar()
    figures.append(fig)

    # Figure 2: Raw frame, orders masked
    fig = plt.figure(figsize=figsize)
    plt.title('Raw frame, orders masked')
    plt.imshow(data_masked, origin='lower', vmin=0, vmax=np.nanmax(data_masked))
    plt.colorbar()
    figures.append(fig)

    # Figure 3: Background model (first iteration)
    fig = plt.figure(figsize=figsize)
    plt.title('Background model')
    plt.imshow(bkg1_background, origin='lower', vmin=0, vmax=np.nanmax(data_masked))
    plt.colorbar()
    figures.append(fig)

    # Figure 4: Background subtracted data (first iteration)
    fig = plt.figure(figsize=figsize)
    plt.title('Background subtracted data')
    plt.imshow(data_original - bkg1_background, origin='lower', vmin=-30, vmax=30)
    plt.colorbar()
    figures.append(fig)

    # Figure 5: Number of masked pixels in background mesh
    fig = plt.figure(figsize=figsize)
    plt.title('Number of masked pixels in background mesh')
    plt.imshow(bkg1_mesh_nmasked, origin='lower', vmin=0, vmax=box_size**2)
    plt.colorbar()
    figures.append(fig)

    # Figure 6: Background model, correction step
    fig = plt.figure(figsize=figsize)
    plt.title('Background model, correction step')
    plt.imshow(bkg2_background, origin='lower', vmin=-10, vmax=10)
    plt.colorbar()
    figures.append(fig)

    # Figure 7: Final background subtracted data
    fig = plt.figure(figsize=figsize)
    plt.title('Final background subtracted data')
    plt.imshow(result, origin='lower', vmin=-30, vmax=30)
    plt.colorbar()
    figures.append(fig)

    return figures


def plot_residuals(residuals, wavelengths, orders_norm, x_norm,
                   weights=None, plottitle=''):
    """
    Create diagnostic plots for wavelength solution residuals.

    Parameters
    ----------
    residuals : ndarray
        Wavelength residuals in nm.
    wavelengths : ndarray
        Wavelength values in nm.
    orders_norm : ndarray
        Normalised order numbers used for colour mapping.
    x_norm : ndarray
        Normalised X positions along each order.
    weights : ndarray, optional
        Per-point weights for the weighted order mean. Defaults to uniform
        weights if not provided.
    plottitle : str, optional
        Title string for all three subplots.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with 3 subplots: residuals vs wavelength, residual rms per
        order, and residuals vs normalised X position.
    """
    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(bottom=.07, left=0.14, right=0.96, top=0.95, hspace=0.40)
    ax1 = fig.add_subplot(311)
    ax3 = fig.add_subplot(312)
    ax2 = fig.add_subplot(313)

    ax1.set_title(plottitle)
    ax1.set_ylabel('Residuals [m/s]')
    ax1.set_xlabel('Wavelength [nm]')
    ax1.set_ylim(-22, 22)

    ax1.scatter(wavelengths, residuals / wavelengths * 3e8,
                c=orders_norm, marker='o', cmap='nipy_spectral', s=2, rasterized=True)

    ordermeans = []
    orderaverages = []
    orderrmss = []
    wavemeans = []
    unique_orders = np.unique(orders_norm)
    wavelengths = np.asarray(wavelengths)

    if weights is None:
        weights = residuals * 0 + 1.

    for o in unique_orders:
        mask = np.where(orders_norm == o)
        res_ms = (residuals / wavelengths * 3e8)
        oa = np.average(res_ms[mask], weights=weights[mask])
        om = np.nanmean(res_ms[mask])
        wm = np.nanmean(wavelengths[mask])
        os = np.nanstd(res_ms[mask])
        ordermeans.append(om)
        orderaverages.append(oa)
        orderrmss.append(os)
        wavemeans.append(wm)

    res_ms_all = residuals / wavelengths * 3e8
    scatter = np.std(res_ms_all)
    ax1.scatter(wavemeans, ordermeans,    marker='p', c='k', s=50, label='Unweighted mean')
    ax1.scatter(wavemeans, orderaverages, marker='*', c='k', s=50, label='Weighted mean')
    ax1.legend(loc=1)
    ax1.hlines(np.mean(res_ms_all), np.min(wavelengths), np.max(wavelengths))
    ax1.text(0.05, 0.05,
             f'std: {scatter:.1f} m/s, mean: {np.mean(res_ms_all):.2f} m/s',
             transform=ax1.transAxes)

    ax3.set_title(plottitle)
    ax3.set_ylabel('Residual rms [m/s]')
    ax3.set_xlabel('Wavelength [nm]')
    ax3.set_ylim(0, 15)
    ax3.scatter(wavemeans, orderrmss, marker='o', c=unique_orders, s=50,
                cmap='nipy_spectral', rasterized=True)

    ax2.set_title(plottitle)
    ax2.set_ylabel('Residuals [m/s]')
    ax2.set_xlabel('X (normalized)')
    ax2.set_ylim(-22, 22)
    if len(x_norm) == len(residuals):
        ax2.scatter(x_norm, res_ms_all, s=2, c=orders_norm, marker='o',
                    cmap='nipy_spectral', rasterized=True)

    return fig


def plot_fiber_combination(order, wave1, wave2, wave3,
                           intensity1_2, intensity2, intensity3_2,
                           weights1_2, weights2, weights3_2,
                           intensity, error,
                           error1_2, error2, error3_2,
                           median_intensity, kappa_sigma):
    """
    Create diagnostic plot for fiber combination for one order.

    Parameters
    ----------
    order : int
        Order number for title.
    wave1, wave2, wave3 : ndarray
        Wavelength arrays for fibers 2, 3, 4.
    intensity1_2, intensity2, intensity3_2 : ndarray
        Intensity arrays (fibers 2 and 4 interpolated onto the fiber 3 grid).
    weights1_2, weights2, weights3_2 : ndarray
        Weight arrays (1/variance) for each fiber.
    intensity : ndarray
        Combined intensity.
    error : ndarray
        Combined error.
    error1_2, error2, error3_2 : ndarray
        Error arrays for individual fibers.
    median_intensity : ndarray
        Median intensity across fibers.
    kappa_sigma : float
        Sigma clipping threshold used.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with 4 subplots: intensity, sigma deviations, weights, SNR.
    """
    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(bottom=0.07, left=0.14, right=0.96, top=0.95, hspace=0.40)
    ax0 = fig.add_subplot(411)
    ax1 = fig.add_subplot(412, sharex=ax0)
    ax2 = fig.add_subplot(413, sharex=ax0)
    ax3 = fig.add_subplot(414, sharex=ax0)

    ax0.set_title(f'Order: {order}')
    ax3.set_xlabel('Wavelength (nm)')

    # ax0 - Intensity
    ax0.set_ylabel('Intensity (DN)')
    ax0.plot(wave2, intensity1_2, 'r', label='Fiber 2', rasterized=True)
    ax0.plot(wave2, intensity2,   'g', label='Fiber 3', rasterized=True)
    ax0.plot(wave2, intensity3_2, 'b', label='Fiber 4', rasterized=True)
    ax0.plot(wave2, intensity,    'k', label='Combined', rasterized=True)
    ax0.legend()

    # ax1 - Sigma deviations from median
    ax1.set_title('Sigma deviations from median')
    ax1.set_ylabel('Sigma')
    ax1.set_ylim(0, kappa_sigma + 1)
    ax1.plot([np.min(wave2), np.max(wave2)], [kappa_sigma, kappa_sigma], 'k:', rasterized=True)
    ax1.plot(wave2, np.abs(intensity1_2 - median_intensity) / np.sqrt(error1_2), 'r', rasterized=True)
    ax1.plot(wave2, np.abs(intensity2   - median_intensity) / np.sqrt(error2),   'g', rasterized=True)
    ax1.plot(wave2, np.abs(intensity3_2 - median_intensity) / np.sqrt(error3_2), 'b', rasterized=True)

    # ax2 - Weights
    ax2.set_title('Weights')
    ax2.set_ylabel('Weights (1/variance)')
    ax2.plot(wave1, weights1_2, 'r', rasterized=True)
    ax2.plot(wave2, weights2,   'g', rasterized=True)
    ax2.plot(wave3, weights3_2, 'b', rasterized=True)

    # ax3 - SNR
    ax3.set_title('SNR')
    ax3.set_ylabel('SNR')
    ax3.plot(wave2, intensity / np.sqrt(error), 'k', rasterized=True)

    return fig


def plot_calibfiber_offset(xs, shifts, orders, wavelengths, splfits,
                           fig=None, plottitle=''):
    """
    Create diagnostic plots for calibration fiber offset measurements.

    Can be called once to create a new figure, or repeatedly with the same
    ``fig`` to overlay multiple datasets before saving.

    Parameters
    ----------
    xs : ndarray
        Normalised X positions along each order.
    shifts : ndarray
        Measured fiber offsets in m/s.
    orders : ndarray
        Order numbers, used for colour mapping.
    wavelengths : ndarray
        Wavelength values in nm.
    splfits : ndarray
        Spline fit values to the offsets.
    fig : matplotlib.figure.Figure, optional
        Existing figure to append data to. If ``None`` a new figure is
        created with the required axes layout.
    plottitle : str, optional
        Suffix appended to each subplot title.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with 3 subplots: offset vs wavelength with spline fit,
        per-order fit residuals, and offset vs normalised X.
    """
    if fig is None:
        fig = plt.figure(figsize=(8, 10))
        fig.subplots_adjust(bottom=.07, left=0.14, right=0.96, top=0.95, hspace=0.40)
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312, sharex=ax1)
        ax3 = fig.add_subplot(313)

        ax1.set_title('Reference fiber offset ' + plottitle)
        ax1.set_xlabel('Wavelength (nm)')
        ax1.set_ylabel('Offset (m/s)')
        ax1.set_ylim(np.nanmean(shifts) - 30, np.nanmean(shifts) + 30)

        ax2.set_title('Reference fiber offset ' + plottitle)
        ax2.set_xlabel('Wavelength (nm)')
        ax2.set_ylabel('Offset (m/s)')
        ax2.set_ylim(np.nanmin(splfits) - 2, np.nanmax(splfits) + 2)

        ax3.set_title('Reference fiber offset ' + plottitle)
        ax3.set_xlabel('X (normalized)')
        ax3.set_ylabel('Offset (m/s)')
        ax3.set_ylim(np.nanmean(shifts) - 30, np.nanmean(shifts) + 30)

    ax = fig.axes

    ax[0].scatter(wavelengths, shifts, c=orders, cmap='nipy_spectral',
                  rasterized=True, marker='.', s=2)
    ax[2].scatter(xs, shifts, c=orders, cmap='nipy_spectral',
                  rasterized=True, marker='.', s=2)

    wavelenths_o = []
    splfits_mean_o = []
    splfits_std_o = []
    unique_orders = np.unique(orders)

    for o in unique_orders:
        idx = np.where(orders == o)
        ax[0].plot(wavelengths[idx], splfits[idx], c='k', rasterized=True)
        ax[2].plot(xs[idx], splfits[idx], c='k', rasterized=True)
        wavelenths_o  = np.append(wavelenths_o,  np.nanmean(wavelengths[idx]))
        splfits_mean_o = np.append(splfits_mean_o, np.nanmean(splfits[idx]))
        splfits_std_o  = np.append(splfits_std_o,  np.nanstd(splfits[idx]))

    ax[1].set_ylim(np.nanmin(splfits_mean_o) - 1, np.nanmax(splfits_mean_o) + 1)
    ax[1].hlines(0, np.min(wavelengths), np.max(wavelengths))
    ax[1].scatter(wavelenths_o, splfits_mean_o, c=unique_orders,
                  cmap='nipy_spectral', rasterized=True, marker='o')
    ax[1].text(0.05, 0.05,
               f'mean: {np.mean(splfits_mean_o):.2f} m/s', transform=ax[1].transAxes)

    return fig


def plot_etalon_residuals(wavelengths, residuals, orders,
                          plottitle='', plotnumber=0, fig=None):
    """
    Create diagnostic plots for etalon wavelength solution residuals.

    Can be called up to four times with the same ``fig`` (using
    ``plotnumber`` 0–3) to fill all subplots before saving.

    Parameters
    ----------
    wavelengths : ndarray
        Wavelength values in nm.
    residuals : ndarray
        Residuals in m/s.
    orders : ndarray
        Order numbers, used for colour mapping.
    plottitle : str, optional
        Suffix appended to the subplot title.
    plotnumber : int, optional
        Index (0–3) of the subplot to populate. Defaults to 0.
    fig : matplotlib.figure.Figure, optional
        Existing figure to append data to. If ``None`` a new figure is
        created with four shared-x subplots.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with 4 subplots.
    """
    if fig is None:
        fig = plt.figure(figsize=(8, 10))
        fig.subplots_adjust(bottom=.07, left=0.14, right=0.96, top=0.95, hspace=0.40)
        ax1 = fig.add_subplot(411)
        ax2 = fig.add_subplot(412, sharex=ax1)
        ax3 = fig.add_subplot(413, sharex=ax1)
        ax4 = fig.add_subplot(414, sharex=ax1)

    ax = fig.axes

    ax[plotnumber].set_title('Etalon residuals ' + plottitle)
    ax[plotnumber].set_xlabel('Wavelength (nm)')
    ax[plotnumber].set_ylabel('Residuals (m/s)')
    ax[plotnumber].scatter(wavelengths, residuals, c=orders, cmap='nipy_spectral',
                           rasterized=True, marker='.', s=2)
    ax[plotnumber].set_ylim(np.nanmean(residuals) - 30, np.nanmean(residuals) + 30)
    ax[plotnumber].text(0.6, 0.05,
                        f'std: {np.nanstd(residuals):.1f} m/s, mean: {np.nanmean(residuals):.2f} m/s',
                        transform=ax[plotnumber].transAxes)

    return fig
