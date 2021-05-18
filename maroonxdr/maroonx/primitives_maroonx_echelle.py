#
#                                                                       DRAGONS
#
#                                                 primitives_maroonx_echelle.py
# ------------------------------------------------------------------------------

from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.gemini import gemini_tools as gt
import numpy as np
from scipy import ndimage

from geminidr.core import Spect
from .primitives_maroonx import MAROONX
from . import parameters_maroonx_echelle

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

@parameter_override
class MAROONXEchelle(MAROONX, Spect):
    """
    This class contains primitives that applies to all MAROON-X echelle
    data.
    """

    tagset = {'GEMINI', 'MAROONX', 'ECHELLE', 'SPECT'}

    def __init__(self, adinputs, **kwargs):
        super(MAROONXEchelle, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'maroonxdr.maroonx.lookups'
        self._param_update(parameters_maroonx_echelle)

    def myNewPrimitive(self, adinputs=None, **params):
        """
        Description...
        
        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        param2: blah
            blah, blah

        Returns
        -------
        """

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Get params out
        param2 = params['param2']

        # Initialize the list of output AstroData objects
        # It is also possible to modify adinputs in place.
        adoutputs = []

        for ad in adinputs:

            # Do whatever checks on the input are necessary, for example:
            # Check whether this primitive as been run already.
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has"
                            "already been processed by myNewPrimitive".
                            format(ad.filename))
                continue

            # -----------------------
            # DR algorithm goes here
            # -----------------------

            # Timestamp
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

            adoutputs.append(ad_out)

        return adoutputs

    def correctImageOrientation(self, adinputs=None, debug_level=0):
        """
        Correct image orientation to proper echelle format.

        flips image so that left lower corner is bluest wavelength, upper right corner is reddest wavelength.
        Echelle orders go from left to right.

        Args:
            adinputs --- previously --- img (np.ndarray): input 2d echelle spectrum
            debug_level (int): debug level
        Returns:
            np.ndarray: image with correct orientation
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:

            if debug_level > 0:
                plt.figure()
                plt.title('Original Image')
                plt.imshow(ad.data[0], origin='lower')

            if ad.image_orientation()[0]:  # flip up-down
                ad.data[0] = np.flipud(ad.data[0])
            if ad.image_orientation()[1]:  # flip left-right (dispersion direction)
                ad.data[0] = np.fliplr(ad.data[0])

            if debug_level > 0:
                plt.figure()
                plt.title('Orientation Corrected image')
                plt.imshow(ad.data[0], origin='lower')
                plt.show()

        # return imgs
        gt.mark_history(adinputs, primname=self.myself(), keyword=timestamp_key)

        return adinputs


    def findStripes(self, adinput=None, badpixelmap=None, deg_polynomial=5, median_filter=1, gauss_filter_sigma=3., min_peak=0.25,
                     debug_level=0):
        """
        Locates and fits stripes in a flat field echelle spectrum.

        Starting in the central column, the algorithm identifies peaks and traces each stripe to the edge of the detector
        by following the brightest pixels along each order. It then fits a polynomial to each stripe. To improve
        algorithm stability, the image is first median filtered and then smoothed with a gaussian. It not only eliminates
        noise, but also ensures that the cross section profile of the flat becomes peaked in the middle, which helps to
        identify the center of each stripe. Choose gauss_filter accordingly. To avoid false positives, only peaks above a
        certain (relative) intensity threshold are used.

        Args:
            adinput -- previously -- flat (np.ndarray): dark corrected flat field spectrum
            badpixelmap (np.ndarray): bad pixel map (optional)
            deg_polynomial (int): degree of the polynomial fit
            median_filter (int): median filter parameter
            gauss_filter_sigma (float): sigma of the gaussian filter used to smooth the image.
            min_peak (float): minimum relative peak height
            debug_level (int): debug level

        Returns:
            list: list of polynomial fits (np.poly1d)
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        ny, nx = adinput.data[0].shape

        # smooth image slightly for noise reduction
        filtered_flat = ndimage.median_filter(adinput.data[0], median_filter)
        filtered_flat = ndimage.gaussian_filter(filtered_flat, gauss_filter_sigma)

        # find peaks in center column
        data = filtered_flat[:, int(nx / 2)]
        peaks = np.r_[True, data[1:] >= data[:-1]] & np.r_[data[:-1] > data[1:], True]

        idx = np.logical_and(peaks, data > min_peak * np.max(data))
        maxima = np.arange(ny)[idx]

        # filter out maxima too close to the boundary to avoid problems
        maxima = maxima[maxima > 3]
        maxima = maxima[maxima < ny - 3]

        if debug_level > 2:
            plt.figure()
            plt.title('Local maxima')
            plt.plot(data)
            plt.plot(np.arange(ny)[peaks], data[peaks], 'r+')
            plt.plot(maxima, data[maxima], 'go')
            plt.show()

        n_order = len(maxima)
        log.info(f'Number of stripes found: {n_order}')

        orders = np.zeros((n_order, nx))
        # walk through to the left and right along the maximum of the order
        for m, row in enumerate(maxima):
            column = int(nx / 2)
            orders[m, column] = row
            start_row = row
            # walk right
            while column + 1 < nx:
                column += 1
                args = [start_row]
                if start_row - 1 > 1:
                    args.append(start_row - 1)
                else:
                    args.append(1)
                if start_row + 1 < ny - 1:
                    args.append(start_row + 1)
                else:
                    args.append(ny - 1)
                p = filtered_flat[args, column]
                # new maximum
                start_row = args[int(np.argmax(p))]
                orders[m, column] = start_row

            column = int(nx / 2)
            start_row = row
            # walk left
            while column > 0:
                column -= 1
                args = [start_row]
                if start_row - 1 > 1:
                    args.append(start_row - 1)
                else:
                    args.append(1)
                if start_row + 1 < ny - 1:
                    args.append(start_row + 1)
                else:
                    args.append(ny - 1)

                p = filtered_flat[args, column]
                # new maximum
                start_row = args[int(np.argmax(p))]
                orders[m, column] = start_row
        # do Polynomial fit for each order
        log.info(f'Fit polynomial of order {deg_polynomial} to each stripe')
        xx = np.arange(nx)
        polynomials = [np.poly1d(np.polyfit(xx, o, deg_polynomial)) for o in orders]

        if debug_level > 0:
            plt.figure()
            plt.title('Traced stripes')
            plt.imshow(filtered_flat, interpolation='none', vmin=np.min(adinput.data[0]), vmax=0.5 * np.max(adinput.data[0]),
                       cmap=plt.get_cmap('gray'), origin='lower')
            for p in polynomials:
                plt.plot(xx, p(xx), 'g', alpha=0.5)
            plt.ylim((0, ny))
            plt.xlim((0, nx))
            plt.show()

        return polynomials