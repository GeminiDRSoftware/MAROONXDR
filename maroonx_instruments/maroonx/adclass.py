import re
from astrodata import (
    Section,
    TagSet,
    astro_data_descriptor,
    astro_data_tag,
    returns_list,
)
from astropy.time import Time, TimeDelta
from gemini_instruments.gemini import AstroDataGemini

from . import lookup

# gemini_keyword_names = dict(overscan_section = 'BIASSEC')

# Define the fiber types that are used to identify the type of data
# These names and fiber configuration are taken from the file:
# 'MAROON-X Data Archiving Notes.pdf'
DARK = 'Dark'
FLAT = 'Flat lamp'
SKY = 'Sky'
OBJECT = 'Target'
ETALON = 'Etalon'
THAR = 'ThAr'
LFC = 'LFC'
IODINE = 'Iodine cell'

FLAT_FIBER_SETUPS = [
    [FLAT, DARK, DARK, DARK, FLAT],
    [DARK, FLAT, FLAT, FLAT, DARK],
    [DARK, DARK, DARK, DARK, FLAT],
]

DARK_FIBER_SETUPS = [
    [DARK, DARK, DARK, DARK, ETALON],
]

SCIENCE_FIBER_SETUPS = [
    [SKY, OBJECT, OBJECT, OBJECT, ETALON],
]

THAR_FIBER_SETUPS = [
    [DARK, THAR, THAR, THAR, THAR],
    [DARK, THAR, THAR, THAR, ETALON],
    [DARK, THAR, THAR, THAR, DARK],
]

ETALON_FIBER_SETUPS = [
    [DARK, ETALON, ETALON, ETALON, ETALON],
    [DARK, ETALON, ETALON, ETALON, IODINE],
]

LFC_FIBER_SETUPS = [
    [DARK, LFC, LFC, LFC, LFC],
    [DARK, LFC, LFC, LFC, ETALON],
    [DARK, LFC, LFC, LFC, DARK],
    [DARK, ETALON, ETALON, ETALON, LFC],
]

WAVECAL_FIBER_SETUPS = THAR_FIBER_SETUPS + ETALON_FIBER_SETUPS + LFC_FIBER_SETUPS


class AstroDataMAROONX(AstroDataGemini):
    # single keyword mapping.  add only the ones that are different
    # from what's already defined in AstroDataGemini.

    __keyword_dict = dict()

    @staticmethod
    def _matches_data(source):
        return (
            source[0].header.get('INSTRUME', '').upper() == 'MAROON-X'
        )  # TODO: add to all headers

    # def _matches_data(source):
    #     if 'HIERARCH MAROONX PUPILCAMERA STATUS' in source[0].header:
    #         return True
    # ---------------
    # Tag definitions
    # ---------------
    @astro_data_tag
    def _tag_instrument(self):
        return TagSet(['MAROONX'])

    # @astro_data_tag
    # def _tag_spect(self):
    #     return TagSet(['SPECT'])

    # @astro_data_tag
    # def _tag_echelle(self):
    #     return TagSet(['ECHELLE'])

    @astro_data_tag
    def _tag_arm(self):
        """Tag the data as either BLUE, RED, or BUNDLE."""
        if self.is_single:
            if self.hdr.get('ARM') == 'BLUE':
                return TagSet(['BLUE'])
            if self.hdr.get('ARM') == 'RED':
                return TagSet(['RED'])

        elif len(self.indices) == 2:
            if self[0].hdr.get('ARM') == 'BLUE' and self[1].hdr.get('ARM') == 'RED':
                return TagSet(['BUNDLE'])
        elif len(self.indices) == 1:
            if not self.filename[0].isdigit():
                # some files have a single arm extension and should
                # be tagged as BUNDLE for the debundle recipe
                return TagSet(['BUNDLE'])
            
            if self[0].hdr.get('ARM') == 'BLUE':
                return TagSet(['BLUE'])
            if self[0].hdr.get('ARM') == 'RED':
                return TagSet(['RED'])
        else:
            return TagSet(['UNDEFINED'])

    @astro_data_tag
    def _tag_exptime(self):
        if self.is_single:
            return TagSet([f'{int(self.hdr.get("EXPTIME"))}s'])
        if len(self.indices) == 1:
            return TagSet([f'{int(self[0].hdr.get("EXPTIME"))}s'])

    @astro_data_tag
    def _tag_dark(self):
        if self.fiber_setup() in DARK_FIBER_SETUPS:
            if hasattr(self[0], 'COEFF_Z0'):
                return TagSet(['DARK', 'DARK_COEFF', 'CAL'])
            if self.phu.get('OBSTYPE') == 'OBJECT':
                return TagSet(['DARK', 'DARK_SYNTH', 'CAL'])
            return TagSet(['DARK', 'CAL'])

    @astro_data_tag
    def _tag_flat(self):
        if self.fiber_setup() in FLAT_FIBER_SETUPS:
            return TagSet(['FLAT', 'CAL'])

    @astro_data_tag
    def _tag_science(self):
        if self.fiber_setup() in SCIENCE_FIBER_SETUPS:
            return TagSet(['SCI', 'SPECT'])

    # @astro_data_tag
    # def _tag_wavecal(self):
    #     # if self.phu.get('FIBER1') == 'Etalon' or self.phu.get('FIBER2') == 'Etalon':
    #     #     if self.phu.get('FIBER5') == 'Etalon':
    #     #         return TagSet(['WAVECAL', 'SPECT', 'CAL'])
    #     if self.fiber_setup() in WAVECAL_FIBER_SETUPS:
    #         return TagSet(['WAVECAL', 'SPECT', 'CAL'])

    @astro_data_tag
    def _tag_etalon(self):
        if self.fiber_setup() in ETALON_FIBER_SETUPS:
            return TagSet(['WAVECAL', 'SPECT', 'ETALON', 'CAL'])        

    @astro_data_tag
    def _tag_thar(self):
        # if (
        #     self.phu.get('FIBER1') == 'ThAr' or self.phu.get('FIBER2') == 'ThAr'
        # ) and self.phu.get('FIBER5') == 'ThAr':
        #     return TagSet(['WAVECAL', 'SPECT', 'ThAr', 'CAL'])
        if self.fiber_setup() in THAR_FIBER_SETUPS:
            return TagSet(['WAVECAL', 'SPECT', 'ThAr', 'CAL'])

    @astro_data_tag
    def _tag_lfc(self):
        # if (
        #     self.phu.get('FIBER1') == 'LFC' or self.phu.get('FIBER2') == 'LFC'
        # ) and self.phu.get('FIBER5') == 'LFC':
        #     return TagSet(['WAVECAL', 'SPECT', 'LFC', 'CAL'])
        if self.fiber_setup() in LFC_FIBER_SETUPS:
            return TagSet(['WAVECAL', 'SPECT', 'LFC', 'CAL'])

    @astro_data_tag
    def _tag_bpm(self):
        if self.phu.get('OBSTYPE') == 'BPM':
            return TagSet(['BPM'])

    # Adapt as required.
    # More tags needs to be added by the MAROON-X DR team
    @astro_data_tag
    def _status_processed_maroonx_cals(self):
        """
        Define the 'processed data' tag set for MAROON-X data.
        """
        kwords = {'PRWAVECAL'}
        if set(self.phu) & kwords:
            return TagSet(['PROCESSED'])

    # ----------------------
    # Descriptor definitions
    # ----------------------

    @astro_data_descriptor
    def instrument(self, generic=False):
        """
        Remove the "-" in the name so that it matches the directories in
        maroonxdr.

        Returns
        -------
        str

        """
        return super().instrument().replace('-', '')

    @returns_list
    @astro_data_descriptor
    def array_name(self):
        """
        Returns a list of the names of the arrays of the extensions, or
        a string if called on a single-extension slice

        Returns
        -------
        list/str
            names of the arrays
        """
        if 'BLUE' in self.tags:
            arrays = [lookup.array_name_b]
        if 'RED' in self.tags:
            arrays = [lookup.array_name_r]
        if 'BUNDLE' in self.tags:
            arrays = [lookup.array_name_b, lookup.array_name_r]

        return arrays

    @astro_data_descriptor
    def fiber_setup(self, short=False):
        """
        Returns the 5 fiber setup for the observation as a list of str

        Parameters
        ----------
        short : bool
            If True, returns short pattern extracted from filename
            (e.g., 'DFFFD' from '20241114T181028Z_DFFFD_b_0008.fits')

        Returns
        -------
        list of str or str
            If short=False: list of fiber names
            If short=True: pattern string extracted from filename
        """
        fibers = [
            self.phu.get('FIBER1'),
            self.phu.get('FIBER2'),
            self.phu.get('FIBER3'),
            self.phu.get('FIBER4'),
            self.phu.get('FIBER5'),
        ]

        if short:
            # Extract pattern from filename
            # Expected format: 
            #  YYYYMMDDTHHMMSSZ_PATTERN_arm_exptime.fits
            # Match pattern between timestamp and arm (b/r)
            filename = self.filename
            match = re.search(r'_([A-Z]{5})_[br]_', filename)
            if match:
                return match.group(1)
            else:
                # Fallback to first letters
                return ''.join([f[0].upper() if f else 'X' for f in fibers])
        return fibers

    @astro_data_descriptor
    def overscan_section(self, pretty=False):
        """
        Returns the overscan (or bias) section.

        Returns
        -------
        list of stings
            Position of the overscan sections using 0-based coordinates.
        """
        ampname = self.array_name()

        if pretty:
            if self.is_single:
                return [lookup.bias_section[amp] for amp in ampname]
            allext = []
            for extampname in ampname:
                allext.append([lookup.bias_section[amp] for amp in extampname])
            return allext

        if self.is_single:
            return [Section.from_string(lookup.bias_section[amp]) for amp in ampname]

        allext = []
        for extampname in ampname:
            allext.append(
                [Section.from_string(lookup.bias_section[amp]) for amp in extampname]
            )
        return allext

    @astro_data_descriptor
    def subtract_overscan_section(self, pretty=False):
        """
        Returns the overscan (or bias) section used for overscan subtraction.

        Returns
        -------
        list of stings
            Position of the overscan sections using 0-based coordinates.
        """
        bs_section = lookup.bias_subtraction_section
        
        # amp names are the keys of the dictionary
        ampname = [list(bs_section.keys()) for _ in self.indices]
        
        if pretty:
            if self.is_single:
                return [bs_section[amp] for amp in ampname[0]]
            allext = []
            for extampname in ampname:
                allext.append([bs_section[amp] for amp in extampname])
            return allext

        if self.is_single:
            return [Section.from_string(bs_section[amp]) for amp in ampname[0]]

        allext = []
        for extampname in ampname:
            allext.append(
                [Section.from_string(bs_section[amp]) for amp in extampname]
            )
        return allext

    @astro_data_descriptor
    def data_section(self, pretty=False):
        """
        Returns the sky-exposable data pixels (or bias) section.

        Returns
        -------
        list of stings
            Position of the sky-exposable sections using 0-based coordinates.
        """
        ampname = self.array_name()

        if pretty:
            if self.is_single:
                return [lookup.data_section[amp] for amp in ampname]
            allext = []
            for extampname in ampname:
                allext.append([lookup.data_section[amp] for amp in extampname])
            return allext
        if self.is_single:
            return [Section.from_string(lookup.data_section[amp]) for amp in ampname]
        allext = []
        for extampname in ampname:
            allext.append(
                [Section.from_string(lookup.data_section[amp]) for amp in extampname]
            )
        return allext

    @astro_data_descriptor
    def array_section(self, pretty=False):
        """
        Returns the array (full amplifier including overscan) sections.  A
        list of strings of 0-based coordinates is returned.

        Returns
        -------
            list of stings
            Position of the array sections using 0-based coordinates.
        """
        ampname = self.array_name()

        if pretty:
            if self.is_single:
                return [lookup.array_section[amp] for amp in ampname]
            allext = []
            for extampname in ampname:
                allext.append([lookup.array_section[amp] for amp in extampname])
            return allext
        if self.is_single:
            return [Section.from_string(lookup.array_section[amp]) for amp in ampname]
        allext = []
        for extampname in ampname:
            allext.append(
                [Section.from_string(lookup.array_section[amp]) for amp in extampname]
            )
        return allext

    @astro_data_descriptor
    def array_subtract_overscan_section(self, pretty=False):
        """
        Returns the array sections that will be corrected by overscan.  A
        list of strings of 0-based coordinates is returned.

        Returns
        -------
            list of stings
            Position of the array sections using 0-based coordinates.
        """
        as_section = lookup.array_subtraction_section
        
        # amp names are the keys of the dictionary
        ampname = [list(as_section.keys()) for _ in self.indices]

        if pretty:
            if self.is_single:
                return [as_section[amp] for amp in ampname[0]]
            allext = []
            for extampname in ampname:
                allext.append([as_section[amp] for amp in extampname])
            return allext
        if self.is_single:
            return [Section.from_string(as_section[amp]) for amp in ampname[0]]
        allext = []
        for extampname in ampname:
            allext.append(
                [Section.from_string(as_section[amp]) for amp in extampname]
            )
        return allext


    @astro_data_descriptor
    def detector_section(
        self, pretty=False
    ):  # only used in BPM ext call as of 10-28-22
        """
        Returns the full frame covered by the detector(s) (all amplifier including overscans).  A
        list of strings of 0-based coordinates is returned.

        Returns
        -------
            list of stings
            Position of the array sections using 0-based coordinates.
        """
        if pretty:
            if self.is_single:
                return lookup.detector_section['RB']
            return [lookup.detector_section['RB']]
        if self.is_single:
            return Section.from_string(lookup.detector_section['RB'])
        return [Section.from_string(lookup.detector_section['RB'])]

    @astro_data_descriptor
    def read_noise(self):
        # TODO: check if read noise is requested in variance or not
        # the problem is that the lookup file and the original pipeline quote the
        # value in data units and in variance,
        # but the header has the value in electrons and not in variance

        # read from header - header is in e- and NOT the variance
        # def _read_noise_variance(hdr):
        #     # read noise variance in data units
        #     return (hdr.get('RDNOISE') / hdr.get('GAIN')) ** 2
        #
        # if self.is_single:
        #     return [_read_noise_variance(self.hdr)]
        # else:
        #     return [_read_noise_variance(self[ext].hdr) for ext in self.indices]

        # old implementation - read noise is in variance
        ampname = self.array_name()
        if self.is_single:
            return [lookup.read_noise[amp] for amp in ampname]
        allext = []
        for extampname in ampname:
            allext.append([lookup.read_noise[amp] for amp in extampname])
        return allext

    @astro_data_descriptor
    def gain(self):
        # read from header
        # if self.is_single:
        #     return [self.hdr.get('GAIN')]
        # else:
        #     return [self[ext].hdr.get('GAIN') for ext in self.indices]

        # old implementation
        ampname = self.array_name()
        if self.is_single:
            return [lookup.gain[amp] for amp in ampname]
        allext = []
        for extampname in ampname:
            allext.append([lookup.gain[amp] for amp in extampname])
        return allext

    @astro_data_descriptor
    def filter_orientation(
        self,
    ):  # this needs to be checked for all analysis utilizing fifth fiber data
        # i.e. dark creation (some value > 0), flat creation (always 0), science extractions (same as dark)
        nd_pos = self.hdr.get('HIERARCH MAROONX ND POSITION')[0]
        return {'ND': nd_pos}

    @astro_data_descriptor
    def image_orientation(self):  # dictionary descriptor
        return {
            'horizontal orientation flip': self.phu.get(
                'HIERARCH MAROONX IMAGE ORIENTATION HORIZONTAL FLIP'
            ),
            'vertical orientation flip': self.phu.get(
                'HIERARCH MAROONX IMAGE ORIENTATION VERTICAL FLIP'
            ),
        }

    @astro_data_descriptor
    def telescope_mjd(self, pretty=False):
        """
        Returns the MJD of the observation as read from the header.

        Parameters
        ----------
        pretty : bool, optional
            If True, returns a Time object.  Default is False.
        Returns
        -------
        float or Time
            MJD of the observation.
        """

        mjd = self.hdr.get('HIERARCH MAROONX TELESCOPE MJD')
        if not pretty:
            return mjd

        if self.is_single:
            return Time(float(mjd), format='mjd', scale='utc')
        else:
            return [Time(float(mjd[i]), format='mjd', scale='utc') for i in self.indices]

    @astro_data_descriptor
    def exposure_time(self, pretty=False):
        """
        Returns the exposure time in seconds.

        Parameters
        ----------
        pretty : bool, optional
            If True, returns a TimeDelta object.  Default is False.

        Returns
        -------
        float or TimeDelta
            Exposure time.
        """
        exposure_time = self.phu.get(self._keyword_for('exposure_time'), -1)
        if exposure_time < 0:
            return None
        
        if not pretty:
            return exposure_time
        return TimeDelta(exposure_time, format='sec')


    @astro_data_descriptor
    def detector_x_bin(self):
        return 1

    @astro_data_descriptor
    def detector_y_bin(self):
        return 1

    # =================================================================
    # Post processing descriptors
    # =================================================================
    # The following descriptors are used to extract information
    # from the processed data that is populated by specific primitives.

    @astro_data_descriptor
    def fiber_drifts(self):
        """
        Returns the 5 fiber drift values as a list of float.
        Drift values are in m/s.

        Note:
        Currently drifts values are calculated for ETALON frames by
        the fitAndApplyEtalonWls primitive.

        Returns
        -------
        list of float
            The drift values for each fiber.
        """
        return [
            self.hdr.get('DRIFT_FIBER_1')[0],
            self.hdr.get('DRIFT_FIBER_2')[0],
            self.hdr.get('DRIFT_FIBER_3')[0],
            self.hdr.get('DRIFT_FIBER_4')[0],
            self.hdr.get('DRIFT_FIBER_5')[0],
        ]