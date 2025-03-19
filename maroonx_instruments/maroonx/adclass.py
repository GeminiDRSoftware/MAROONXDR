from astrodata import Section, TagSet, astro_data_descriptor, astro_data_tag
from gemini_instruments.gemini import AstroDataGemini

from . import lookup

# gemini_keyword_names = dict(overscan_section = 'BIASSEC')

FLAT_FIBER_SETUPS = [
    ['Flat lamp', 'Dark', 'Dark', 'Dark', 'Flat lamp'],
    ['Dark', 'Flat lamp', 'Flat lamp', 'Flat lamp', 'Dark'],
    ['Dark', 'Dark', 'Dark', 'Dark', 'Flat lamp'],
]


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
        if len(self.indices) == 2:
            if self[0].hdr.get('ARM') == 'BLUE' and self[1].hdr.get('ARM') == 'RED':
                return TagSet(['BUNDLE'])
        elif len(self.indices) == 1:
            if self[0].hdr.get('ARM') == 'BLUE':
                return TagSet(['BLUE'])
            if self[0].hdr.get('ARM') == 'RED':
                return TagSet(['RED'])
        else:
            return TagSet(['UNDEFINED'])

    @astro_data_tag
    def _tag_dark(self):
        if self.phu.get('FIBER1') == 'Dark' and self.phu.get('FIBER2') == 'Dark':
            return TagSet(['DARK', 'CAL'])

    @astro_data_tag
    def _tag_science(self):
        if (
            self.phu.get('FIBER1') == 'Sky'
            and self.phu.get('FIBER2')
            == self.phu.get('FIBER3')
            == self.phu.get('FIBER4')
            == 'OBJECT'
            and self.phu.get('FIBER5') == 'Etalon'
        ):
            return TagSet(['SCI', 'SPECT'])

    @astro_data_tag
    def _tag_flat(self):
        if self.fiber_setup() in FLAT_FIBER_SETUPS:
            return TagSet(['FLAT', 'CAL'])

    @astro_data_tag
    def _tag_wavecal(self):
        if self.phu.get('FIBER1') == 'Etalon' or self.phu.get('FIBER2') == 'Etalon':
            if self.phu.get('FIBER5') == 'Etalon':
                return TagSet(['WAVECAL', 'SPECT', 'CAL'])

    @astro_data_tag
    def _tag_thar(self):
        if (
            self.phu.get('FIBER1') == 'ThAr' or self.phu.get('FIBER2') == 'ThAr'
        ) and self.phu.get('FIBER5') == 'ThAr':
            return TagSet(['WAVECAL', 'SPECT', 'ThAr', 'CAL'])

    @astro_data_tag
    def _tag_lfc(self):
        if (
            self.phu.get('FIBER1') == 'LFC' or self.phu.get('FIBER2') == 'LFC'
        ) and self.phu.get('FIBER5') == 'LFC':
            return TagSet(['WAVECAL', 'SPECT', 'LFC', 'CAL'])

    @astro_data_tag
    def _tag_bpm(self):
        if self.phu.get('OBSTYPE') == 'BPM':
            return TagSet(['BPM'])

    # Adapt as required.
    # More tags needs to be added by the MAROON-X DR team

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
            arrays = lookup.array_name_b
        elif 'RED' in self.tags:
            arrays = lookup.array_name_r

        if self.is_single:
            return arrays
        return [arrays]

    @astro_data_descriptor
    def fiber_setup(self):
        """
        Returns the 5 fiber setup for the observation as a list of str

        Returns
        -------
        list of str

        """
        return [
            self.phu.get('FIBER1'),
            self.phu.get('FIBER2'),
            self.phu.get('FIBER3'),
            self.phu.get('FIBER4'),
            self.phu.get('FIBER5'),
        ]

    @astro_data_descriptor
    def data_label(self):
        """
        Returns the target name as recorded by telops

        Returns
        -------
        str name
        """
        return self.phu.get('HIERARCH MAROONX TELESCOPE TARGETNAME')

    @astro_data_descriptor
    def observation_id(self):
        """
        Returns the observations program ID as recorded by telops

        Returns
        -------
        str pid
        """
        return self.phu.get('HIERARCH MAROONX TELESCOPE PROGRAMID')

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
        ampname = self.array_name()
        if self.is_single:
            return [lookup.read_noise[amp] for amp in ampname]
        allext = []
        for extampname in ampname:
            allext.append([lookup.read_noise[amp] for amp in extampname])
        return allext

    @astro_data_descriptor
    def gain(self):
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
        return {'ND': '{:.1f}'.format(self.phu.get('HIERARCH MAROONX ND POSITION'))}

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
    def detector_x_bin(self):
        return 1

    @astro_data_descriptor
    def detector_y_bin(self):
        return 1

    # For a list of expected descriptors, see the appendix in the Astrodata
    # User Manual.
