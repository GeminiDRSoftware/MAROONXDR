from astrodata import astro_data_tag, astro_data_descriptor, returns_list, TagSet
from gemini_instruments import gmu
from gemini_instruments.common import Section
from . import lookup
from gemini_instruments.gemini import AstroDataGemini

import re

class AstroDataMAROONX(AstroDataGemini):

    # single keyword mapping.  add only the ones that are different
    # from what's already defined in AstroDataGemini.

    __keyword_dict = dict()

    @staticmethod
    def _matches_data(source):
        return source[0].header.get('INSTRUME', '').upper() == 'MAROON-X'  # can we change lookup?
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

    @astro_data_tag  # tag or static method? (never mix blue and red, basically two different instruments)
    def _tag_arm(self):
        if re.findall(r"_\w_\d{4}", self.filename)[0][1] == 'b':
            return TagSet(['BLUE'])
        elif re.findall(r"_\w_\d{4}", self.filename)[0][1] == 'r':
            return TagSet(['RED'])
        else:
            return TagSet([''])

    @astro_data_tag
    def _tag_dark(self):
        if self.phu.get('HIERARCH FIBER1') == self.phu.get('HIERARCH FIBER2') == 'Dark'\
                and self.phu.get('HIERARCH FIBER5') == 'Etalon'\
                and self.phu.get('EXPTIME') >= 600.:  # only dark when all fibers are dark
            return TagSet(['DARK', 'CAL'], blocks=['SPECT', 'ECHELLE'])

    @astro_data_tag
    def _tag_flat(self):
        if self.phu.get('HIERARCH FIBER1') == 'Flat' or self.phu.get('HIERARCH FIBER2') == 'Flat' \
                or self.phu.get('HIERARCH FIBER5') == 'Flat':
            # 2 types in descriptors (science flat and cal fiber flat)
            return TagSet(['FLAT', 'CAL'])

    @astro_data_tag
    def _tag_etalon(self):
        if (self.phu.get('HIERARCH FIBER1') == 'Etalon' or self.phu.get('HIERARCH FIBER2') == 'Etalon' \
                or self.phu.get('HIERARCH FIBER5') == 'Etalon')\
                and self.phu.get('EXPTIME') < 600.:
            return TagSet(['ETALON', 'CAL'])

    @astro_data_tag
    def _tag_thar(self):
        if self.phu.get('HIERARCH FIBER1') == 'ThAr' or self.phu.get('HIERARCH FIBER2') == 'ThAr' \
                or self.phu.get('HIERARCH FIBER5') == 'ThAr':
            return TagSet(['THAR', 'CAL'])

    # @astro_data_tag  # no bias frames
    # def _tag_bias(self):
    #     if self.phu.get('OBSTYPE') == 'BIAS':
    #         return TagSet(['BIAS', 'CAL', 'CCD'], blocks=['IMAGE', 'SPECT'])

    @astro_data_tag
    def _tag_iodine(self):
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
    def overscan_section(self):
        """
        Returns the overscan (or bias) section.

        Returns
        -------
        tuple of integers or list of tuples
            Position of the overscan section using Python slice values.
        string or list of strings
            Position of the overscan section using an IRAF section
            format (1-based).
        """
        sections = lookup.bias_section
        return [sections[quad] for quad in sections]

    @astro_data_descriptor
    def ccd_section(self):
        """
        Returns the entire ccd section.  If pretty is False, a
        tuple of 0-based coordinates is returned with format  .
        If pretty is True, a keyword value is returned without parsing as a
        string.  In this format, the coordinates are generally 1-based.
        One tuple or string is return per extension/array.  If more than one
        array, the tuples/strings are return in a list.  Otherwise, the
        section is returned as a tuple or a string.

        Parameters
        ----------
        pretty : bool
         If True, return the formatted string found in the header.

        Returns
        -------
        tuple of integers or list of tuples
            Position of the ccd sections using Python slice values.
        string or list of strings
            Position of the ccd sections using an IRAF section
            format (1-based).
        """
        sections = lookup.ccd_section
        return [sections[quad] for quad in sections]


    # For a list of expected descriptors, see the appendix in the Astrodata
    # User Manual.

