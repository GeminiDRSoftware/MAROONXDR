from astrodata import astro_data_tag, astro_data_descriptor, returns_list, TagSet
from gemini_instruments import gmu
from gemini_instruments.common import Section, section_to_tuple
from . import lookup
from gemini_instruments.gemini import AstroDataGemini

import re
# gemini_keyword_names = dict(overscan_section = 'BIASSEC')

class AstroDataMAROONX(AstroDataGemini):  # ! will need to overhall when arms are combined to one MEF

    # single keyword mapping.  add only the ones that are different
    # from what's already defined in AstroDataGemini.

    __keyword_dict = dict()

    @staticmethod
    def _matches_data(source):
        return source[0].header.get('INSTRUME', '').upper() == 'MAROON-X'  # !add to all data
    # def _matches_data(source):
    #     if 'HIERARCH MAROONX PUPILCAMERA STATUS' in source[0].header:
    #         return True
    # ---------------
    # Tag definitions
    # ---------------
    @astro_data_tag
    def _tag_instrument(self):
        return TagSet(['MAROONX'])

    @astro_data_tag
    def _tag_spect(self):
        return TagSet(['SPECT'])

    @astro_data_tag
    def _tag_echelle(self):
        return TagSet(['ECHELLE'])

    @astro_data_tag  # !add to headers
    def _tag_arm(self):
        if re.findall(r"_\w_\d{4}", self.filename)[0][1] == 'b':
            return TagSet(['BLUE'])
        elif re.findall(r"_\w_\d{4}", self.filename)[0][1] == 'r':
            return TagSet(['RED'])
        else:
            return TagSet(['UNDEFINED'])

    @astro_data_tag
    def _tag_dark(self):
        if self.phu.get('HIERARCH FIBER1') == self.phu.get('HIERARCH FIBER2') == 'Dark'\
                and self.phu.get('HIERARCH FIBER5') == 'Etalon'\
                and self.phu.get('EXPTIME') >= 60.:  # part of dark correction is for etalon light leak into fiber 4
            return TagSet(['DARK', 'CAL'], blocks=['SPECT', 'ECHELLE'])

    @astro_data_tag
    def _tag_flat(self):
        if self.phu.get('HIERARCH FIBER1') == 'Flat' or self.phu.get('HIERARCH FIBER2') == 'Flat' \
                or self.phu.get('HIERARCH FIBER5') == 'Flat':
            # 2 types based on fiber_setup (science flat and cal fibers flat)
            return TagSet(['FLAT', 'CAL'], blocks=['SPECT', 'ECHELLE'])

    @astro_data_tag
    def _tag_etalon(self):
        if (self.phu.get('HIERARCH FIBER1') == 'Etalon' or self.phu.get('HIERARCH FIBER2') == 'Etalon' \
                or self.phu.get('HIERARCH FIBER5') == 'Etalon')\
                and self.phu.get('EXPTIME') < 60.:  # only different than darks by shorter exposure time
            return TagSet(['ETALON', 'CAL'])

    @astro_data_tag
    def _tag_thar(self): # !combine with iodine as ARC class?
        if self.phu.get('HIERARCH FIBER1') == 'ThAr' or self.phu.get('HIERARCH FIBER2') == 'ThAr' \
                or self.phu.get('HIERARCH FIBER5') == 'ThAr':
            return TagSet(['THAR', 'CAL'])

    # @astro_data_tag  # no bias frames
    # def _tag_bias(self):
    #     if self.phu.get('OBSTYPE') == 'BIAS':
    #         return TagSet(['BIAS', 'CAL', 'CCD'], blocks=['IMAGE', 'SPECT'])

    @astro_data_tag
    def _tag_iodine(self):  # !combine with thar as ARC class?
        if self.phu.get('HIERARCH FIBER5') == 'Iodine':  # just simcal
            return TagSet(['IODINE', 'CAL'])

    # The tags above are examples (except the "instrument" tag).  
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
    def array_name(self):  # ! add info to headers and change here to reflect direct access?
        """
        Returns a list of the names of the arrays of the extensions, or
        a string if called on a single-extension slice
        Returns
        -------
        list/str
            names of the arrays
        """
        if 'BLUE' in self.tags:
            return lookup.array_name_b
        elif 'RED' in self.tags:
            return lookup.array_name_r



    @astro_data_descriptor
    def fiber_setup(self):
        """
        Returns the 5 fiber setup for the observation as a list of str

        Returns
        -------
        list of str

        """
        return [self.phu.get('HIERARCH FIBER1'),self.phu.get('HIERARCH FIBER2'),self.phu.get('HIERARCH FIBER3'),
                self.phu.get('HIERARCH FIBER4'),self.phu.get('HIERARCH FIBER5')]

    @astro_data_descriptor
    def overscan_section(self, pretty=False):  # ! add info to headers and change here to reflect direct access?
        """
        Returns the overscan (or bias) section.

        Returns
        -------
        list of stings
            Position of the overscan sections using 0-based coordinates.
        """
        #return self._parse_section(self._keyword_for('overscan_section'), pretty)
        if pretty:
            return [lookup.bias_section[amp] for amp in self.array_name()]
        else:
            return [section_to_tuple(lookup.bias_section[amp]) for amp in self.array_name()]

    @astro_data_descriptor
    def data_section(self, pretty=False):  # ! add info to headers and change here to reflect direct access?
        """
        Returns the sky-exposable data pixels (or bias) section.

        Returns
        -------
        list of stings
            Position of the sky-exposable sections using 0-based coordinates.
        """
        if pretty:
            return [lookup.data_section[amp] for amp in self.array_name()]
        else:
            return [section_to_tuple(lookup.data_section[amp]) for amp in self.array_name()]

    @astro_data_descriptor
    def array_section(self, pretty):  # ! add info to headers and change here to reflect direct access?
        """
        Returns the array (full amplifier including overscan) sections.  A
        list of strings of 0-based coordinates is returned.

        Returns
        -------
            list of stings
            Position of the array sections using 0-based coordinates.
        """
        if pretty:
            return [lookup.array_section[amp] for amp in self.array_name()]
        else:
            return [section_to_tuple(lookup.array_section[amp]) for amp in self.array_name()]

    @astro_data_descriptor
    def read_noise(self):
        return lookup.read_noise

    @astro_data_descriptor
    def gain(self):
        return lookup.gain

    @astro_data_descriptor
    def nd(self):  # reuse filter name? this needs to be checked for all analysis utilizing fifth fiber data
                    # i.e. dark creation (some value > 0), flat creation (always 0), science extractions (same as dark)
        return '{:.1f}'.format(self.phu.get("HIERARCH MAROONX ND POSITION"))

    # For a list of expected descriptors, see the appendix in the Astrodata
    # User Manual.

