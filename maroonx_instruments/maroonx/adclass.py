"""
AstroData instrument definition for the MAROON-X spectrograph.

This module defines :class:`AstroDataMAROONX`, the AstroData class that the
DRAGONS recipe system uses to recognize and classify MAROON-X frames. Its
``@astro_data_tag`` methods assign the tag set (arm, frame type, processing
status) that selects recipes, and its ``@astro_data_descriptor`` methods
expose header keywords and static detector properties under stable names.

Frame types are identified by the illumination pattern of the five fibers,
read from the ``FIBER1`` to ``FIBER5`` header keywords. The module-level
``*_FIBER_SETUPS`` lists enumerate the recognized patterns for each frame
type; the fiber names follow the document 'MAROON-X Data Archiving Notes'.
Detector geometry, gain, and read noise come from the static tables in
:mod:`~maroonx_instruments.maroonx.lookup`.
"""
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
    """
    AstroData class for frames taken with the MAROON-X spectrograph.

    Files are matched to this class when their primary header carries
    ``INSTRUME = 'MAROON-X'``. The tag methods classify each frame by arm
    (``BLUE``, ``RED``, or ``BUNDLE`` for the two-arm archive files), by
    frame type derived from the five-fiber illumination pattern (``DARK``,
    ``FLAT``, ``SCI``, ``ETALON``, ``ThAr``, ``LFC``), and by processing
    status (``PROCESSED``, ``BARYCOR``). The descriptors expose header
    keywords (fiber setup, exposure time, telescope MJD, ND filter
    position) and per-amplifier detector constants (gain, read noise,
    section geometry) served from the
    :mod:`~maroonx_instruments.maroonx.lookup` tables.
    """

    # single keyword mapping.  add only the ones that are different
    # from what's already defined in AstroDataGemini.

    __keyword_dict = dict()

    @staticmethod
    def _matches_data(source):
        """Match files whose primary header has ``INSTRUME = 'MAROON-X'``."""
        return (
            source[0].header.get('INSTRUME', '').upper() == 'MAROON-X'
        )

    # ---------------
    # Tag definitions
    # ---------------
    @astro_data_tag
    def _tag_instrument(self):
        """Tag all frames as MAROONX."""
        return TagSet(['MAROONX'])

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
        """Tag the frame with its integer exposure time, e.g. ``60s``."""
        if self.is_single:
            return TagSet([f'{int(self.hdr.get("EXPTIME"))}s'])
        if len(self.indices) == 1:
            return TagSet([f'{int(self[0].hdr.get("EXPTIME"))}s'])
        elif len(self.indices) == 2:
            exptimes = set(self.hdr.get("EXPTIME"))
            if len(exptimes) == 1:
                return TagSet([f'{int(exptimes.pop())}s'])

    @astro_data_tag
    def _tag_dark(self):
        """Tag darks as DARK and CAL, with DARK_COEFF or DARK_SYNTH variants."""
        if self.fiber_setup() in DARK_FIBER_SETUPS:
            ext = self if self.is_single else self[0]
            if hasattr(ext, 'COEFF_Z0'):
                return TagSet(['DARK', 'DARK_COEFF', 'CAL'])
            if self.phu.get('OBSTYPE') == 'OBJECT':
                return TagSet(['DARK', 'DARK_SYNTH', 'CAL'])
            return TagSet(['DARK', 'CAL'])

    @astro_data_tag
    def _tag_flat(self):
        """Tag flat frames as FLAT and CAL."""
        if self.fiber_setup() in FLAT_FIBER_SETUPS:
            return TagSet(['FLAT', 'CAL'])

    @astro_data_tag
    def _tag_science(self):
        """Tag science frames as SCI and SPECT."""
        if self.fiber_setup() in SCIENCE_FIBER_SETUPS:
            return TagSet(['SCI', 'SPECT'])

    @astro_data_tag
    def _tag_etalon(self):
        """Tag etalon frames as WAVECAL, SPECT, ETALON, and CAL."""
        if self.fiber_setup() in ETALON_FIBER_SETUPS:
            return TagSet(['WAVECAL', 'SPECT', 'ETALON', 'CAL'])        

    @astro_data_tag
    def _tag_thar(self):
        """Tag ThAr frames as WAVECAL, SPECT, ThAr, and CAL."""
        if self.fiber_setup() in THAR_FIBER_SETUPS:
            return TagSet(['WAVECAL', 'SPECT', 'ThAr', 'CAL'])

    @astro_data_tag
    def _tag_lfc(self):
        """Tag laser frequency comb frames as WAVECAL, SPECT, LFC, and CAL."""
        if self.fiber_setup() in LFC_FIBER_SETUPS:
            return TagSet(['WAVECAL', 'SPECT', 'LFC', 'CAL'])

    @astro_data_tag
    def _tag_bpm(self):
        """Tag bad pixel mask files (``OBSTYPE = 'BPM'``) as BPM."""
        if self.phu.get('OBSTYPE') == 'BPM':
            return TagSet(['BPM'])

    @astro_data_tag
    def _status_processed_maroonx_cals(self):
        """Tag frames with PRWAVECAL or PRDKCOEF as PROCESSED."""
        kwords = {'PRWAVECAL', 'PRDKCOEF'}
        if set(self.phu) & kwords:
            return TagSet(['PROCESSED'])

    @astro_data_tag
    def _tag_barycor(self):
        """Tag barycentric-corrected frames as BARYCOR."""
        if 'BARYCENTRIC_CORRECTION_APPLIED' in self.phu:
            return TagSet(['BARYCOR'])

    # ----------------------
    # Descriptor definitions
    # ----------------------

    @astro_data_descriptor
    def arm(self):
        """
        Return the value of the ``ARM`` keyword for each extension.

        Returns
        -------
        str or list of str
            'BLUE' or 'RED' for a single-extension slice; otherwise one
            value per extension, e.g. ``['BLUE', 'RED']`` for a bundle.
        """
        return self.hdr.get('ARM')

    @astro_data_descriptor
    def camera(self, stripID=False, pretty=False):
        """
        Return the arm name as the camera identifier.

        This descriptor exists so that FitsStorage and the calibration
        database can match calibrations by arm via ``Header.camera``.

        Parameters
        ----------
        stripID : bool
            Inherited from the Gemini descriptor signature; unused.
        pretty : bool
            Inherited from the Gemini descriptor signature; unused.

        Returns
        -------
        str or None
            'BLUE' or 'RED', or None for a bundle with both arms.
        """
        arm = self.hdr.get('ARM')
        if isinstance(arm, list):
            return arm[0] if len(arm) == 1 else None
        return arm

    @astro_data_descriptor
    def instrument(self, generic=False):
        """
        Return the instrument name without the dash, i.e. 'MAROONX'.

        The header value 'MAROON-X' is stripped of its dash so that the
        name matches the ``maroonxdr`` and ``maroonx_instruments``
        package directories used by the recipe system.

        Parameters
        ----------
        generic : bool
            Inherited from the Gemini descriptor signature; unused.

        Returns
        -------
        str
            'MAROONX'
        """
        return super().instrument().replace('-', '')

    @returns_list
    @astro_data_descriptor
    def array_name(self):
        """
        Return the amplifier names of each extension.

        The names are the keys into the
        :mod:`~maroonx_instruments.maroonx.lookup` tables (gain, read
        noise, section geometry): the blue arm reads out through four
        quadrants 'Q1' to 'Q4', the red arm through two halves 'R1' and
        'R2'.

        Returns
        -------
        list
            For a single-extension slice, the list of amplifier names
            (e.g. ``['Q1', 'Q2', 'Q3', 'Q4']``); otherwise one such
            list per extension, with a bundle getting both the blue
            and the red lists.
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
        Return the illumination setup of the five fibers.

        The fiber names are read from the ``FIBER1`` to ``FIBER5``
        keywords of the primary header (e.g. ``['Sky', 'Target',
        'Target', 'Target', 'Etalon']`` for a science frame). The tag
        methods compare this list against the module-level
        ``*_FIBER_SETUPS`` patterns to classify the frame.

        With ``short=True`` the setup is returned as the five-letter
        code used in MAROON-X filenames (one letter per fiber, e.g.
        'SOOOE' for the science setup above), following the naming
        convention of the legacy pipeline.

        Parameters
        ----------
        short : bool
            If True, return the five-letter code parsed from the
            filename (the pattern between the timestamp and the arm
            letter, as in '20241114T181028Z_DFFFD_b_0008.fits'). If
            the filename does not match, fall back to the first
            letters of the header fiber names, with 'X' for a
            missing keyword.

        Returns
        -------
        list of str or str
            List of five fiber names, or the five-letter code if
            ``short=True``.
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
        Return the overscan (bias) section of each amplifier.

        The sections come from the ``bias_section`` table in
        :mod:`~maroonx_instruments.maroonx.lookup`. Note that the
        pipeline's ``subtractOverscan`` primitive does not use this
        descriptor; it works on the legacy quadrant sections from
        :meth:`subtract_overscan_section`.

        Parameters
        ----------
        pretty : bool
            If True, return the sections as the 1-based inclusive
            strings of the lookup table instead of Section objects.

        Returns
        -------
        list of astrodata.Section or list of str
            For a single-extension slice, one section per amplifier;
            otherwise one list per extension.
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
        Return the sections whose mean is used for overscan subtraction.

        These are the legacy overscan regions, one per detector quadrant
        (keys 'RB1' to 'RB4' of the ``bias_subtraction_section`` table
        in :mod:`~maroonx_instruments.maroonx.lookup`). They are the
        same for both arms and reproduce the quadrant arithmetic of the
        legacy pipeline: the ``subtractOverscan`` primitive subtracts
        the mean of each of these regions from the matching quadrant of
        :meth:`array_subtract_overscan_section`.

        Parameters
        ----------
        pretty : bool
            If True, return the sections as the 1-based inclusive
            strings of the lookup table instead of Section objects.

        Returns
        -------
        list of astrodata.Section or list of str
            For a single-extension slice, one section per quadrant;
            otherwise one such list per extension.
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
        Return the exposed-pixel section of each amplifier.

        The sections come from the ``data_section`` table in
        :mod:`~maroonx_instruments.maroonx.lookup` and mark where the
        real (light-sensitive) pixels of each amplifier lie on the raw
        frame, excluding the overscan regions.

        Parameters
        ----------
        pretty : bool
            If True, return the sections as the 1-based inclusive
            strings of the lookup table instead of Section objects.

        Returns
        -------
        list of astrodata.Section or list of str
            For a single-extension slice, one section per amplifier;
            otherwise one list per extension.
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
        Return the section each amplifier occupies after overscan removal.

        The sections come from the ``array_section`` table in
        :mod:`~maroonx_instruments.maroonx.lookup` and give the
        destination of each amplifier's exposed pixels on the trimmed
        frame, once the overscan regions have been removed.

        Parameters
        ----------
        pretty : bool
            If True, return the sections as the 1-based inclusive
            strings of the lookup table instead of Section objects.

        Returns
        -------
        list of astrodata.Section or list of str
            For a single-extension slice, one section per amplifier;
            otherwise one list per extension.
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
        Return the quadrants corrected during overscan subtraction.

        These are the four quadrants of the raw frame (keys 'RB1' to
        'RB4' of the ``array_subtraction_section`` table in
        :mod:`~maroonx_instruments.maroonx.lookup`), the same for both
        arms. The ``subtractOverscan`` primitive subtracts from each
        quadrant the mean of the matching overscan region from
        :meth:`subtract_overscan_section`.

        Parameters
        ----------
        pretty : bool
            If True, return the sections as the 1-based inclusive
            strings of the lookup table instead of Section objects.

        Returns
        -------
        list of astrodata.Section or list of str
            For a single-extension slice, one section per quadrant;
            otherwise one list per extension.
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
    def detector_section(self, pretty=False):
        """
        Return the full raw frame covered by the detector.

        The section is the single '[1:4400, 1:4400]' entry of the
        ``detector_section`` table in
        :mod:`~maroonx_instruments.maroonx.lookup`: the complete raw
        frame including the overscan regions, identical for both arms.
        No MAROON-X code calls this descriptor directly; it is consumed
        by the DRAGONS BPM machinery when ``addDQ`` clips the bad pixel
        mask to the science frame.

        Parameters
        ----------
        pretty : bool
            If True, return the section as the 1-based inclusive
            string of the lookup table instead of a Section object.

        Returns
        -------
        astrodata.Section or str
            The full-frame section for a single-extension slice;
            wrapped in a one-element list otherwise.
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
        """
        Return the read noise of each amplifier, quoted as a variance.

        The values come from the ``read_noise`` table in
        :mod:`~maroonx_instruments.maroonx.lookup` and are the variance
        of the read noise in data units (DN squared), not a noise in
        electrons: 1.14 for every blue quadrant and 1.63 for both red
        halves, matching the values hardcoded in the legacy pipeline.
        ``addVAR`` and the extraction primitives add them directly to
        their variance models. Note that this differs from the
        ``RDNOISE`` header keyword, which is in electrons.

        Returns
        -------
        list of float
            For a single-extension slice, one value per amplifier;
            otherwise one such list per extension.
        """
        ampname = self.array_name()
        if self.is_single:
            return [lookup.read_noise[amp] for amp in ampname]
        allext = []
        for extampname in ampname:
            allext.append([lookup.read_noise[amp] for amp in extampname])
        return allext

    @astro_data_descriptor
    def gain(self):
        """
        Return the gain of each amplifier in electrons per DN.

        The values come from the ``gain`` table in
        :mod:`~maroonx_instruments.maroonx.lookup`: 2.72 for every
        blue quadrant and 2.74 for both red halves, matching the
        values hardcoded in the legacy pipeline.

        Returns
        -------
        list of float
            For a single-extension slice, one value per amplifier;
            otherwise one such list per extension.
        """
        ampname = self.array_name()
        if self.is_single:
            return [lookup.gain[amp] for amp in ampname]
        allext = []
        for extampname in ampname:
            allext.append([lookup.gain[amp] for amp in extampname])
        return allext

    @astro_data_descriptor
    def filter_orientation(self):
        """
        Return the position of the simultaneous-calibration ND filter.

        The value is read from the ``HIERARCH MAROONX ND POSITION``
        keyword of the first extension. The neutral density filter
        attenuates the light fed to the fifth (simultaneous
        calibration) fiber, so frames reduced together must share the
        same setting: darks and science frames are taken with a
        matching, usually nonzero, position and flats with 0. The
        ``checkND`` primitive enforces consistency across an input
        set.

        Returns
        -------
        dict
            One-entry dictionary ``{'ND': position}``.
        """
        nd_pos = self.hdr.get('HIERARCH MAROONX ND POSITION')[0]
        return {'ND': nd_pos}

    @astro_data_descriptor
    def image_orientation(self):
        """
        Return the image flip flags recorded by the instrument.

        The flags are read from the ``HIERARCH MAROONX IMAGE
        ORIENTATION HORIZONTAL FLIP`` and ``VERTICAL FLIP`` keywords
        of the primary header. They record whether the frame needs
        flipping to reach the standard echelle orientation (bluest
        wavelength in the lower left corner). Currently unused: the
        ``correctImageOrientation`` primitive decides the flips from
        the arm tag instead.

        Returns
        -------
        dict
            Keys 'horizontal orientation flip' and 'vertical
            orientation flip' with the header values.
        """
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
        Return the MJD written by the telescope during the exposure.

        The value is read from the ``HIERARCH MAROONX TELESCOPE MJD``
        keyword. The barycentric correction primitive uses it as an
        alternative time base to the UT start time, interpreting it
        (via its ``start_time`` parameter) as written either at the
        start of the exposure or at the end of the readout.

        Parameters
        ----------
        pretty : bool
            If True, return UTC-scale Time objects instead of the
            float header values.

        Returns
        -------
        float or astropy.time.Time
            MJD of the observation; one value (or Time) per extension
            when called on a multi-extension object.
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
        Return the exposure time in seconds.

        The value is read from the ``EXPTIME`` keyword of the primary
        header. None is returned when the keyword is missing or its
        value is negative.

        Parameters
        ----------
        pretty : bool
            If True, return a TimeDelta object instead of the float
            value.

        Returns
        -------
        float or astropy.time.TimeDelta or None
            Exposure time in seconds.
        """
        exposure_time = self.phu.get(self._keyword_for('exposure_time'), -1)
        if exposure_time < 0:
            return None
        
        if not pretty:
            return exposure_time
        return TimeDelta(exposure_time, format='sec')


    @astro_data_descriptor
    def detector_x_bin(self):
        """Return the X-axis binning factor, always 1 (MAROON-X is unbinned)."""
        return 1

    @astro_data_descriptor
    def detector_y_bin(self):
        """Return the Y-axis binning factor, always 1 (MAROON-X is unbinned)."""
        return 1

    # =================================================================
    # Post processing descriptors
    # =================================================================
    # The following descriptors are used to extract information
    # from the processed data that is populated by specific primitives.

    @astro_data_descriptor
    def fiber_drifts(self):
        """
        Return the five per-fiber instrumental drift values in m/s.

        The values are read from the ``DRIFT_FIBER_1`` to
        ``DRIFT_FIBER_5`` keywords of the first extension. They are
        written by the ``fitAndApplyEtalonWLS`` primitive when it
        measures the instrumental drift of an ETALON frame; on frames
        that have not been through that primitive, a list of None
        values is returned.

        Returns
        -------
        list of float
            One drift value per fiber.
        """
        return [
            self.hdr.get('DRIFT_FIBER_1')[0],
            self.hdr.get('DRIFT_FIBER_2')[0],
            self.hdr.get('DRIFT_FIBER_3')[0],
            self.hdr.get('DRIFT_FIBER_4')[0],
            self.hdr.get('DRIFT_FIBER_5')[0],
        ]