"""
Calibration association rules for MAROON-X.

This module defines the CalibrationMAROONX class that FitsStorage uses
to match science frames with their calibrations (flats, darks, BPMs, etc.).

Matching rules (from 'MaroonX calibration and association rules.md'):

- Flat: same arm, closest in time
- Dark: same arm, same exposure time, closest in time
- BPM: same arm, single master
"""

from fits_storage.core.orm.header import Header
from fits_storage.cal.calibration import Calibration


class CalibrationMAROONX(Calibration):
    """
    Calibration matching rules for the MAROON-X spectrograph.

    MAROON-X has no instrument-specific ORM table; all matching
    is done via Header columns (instrument, camera=arm, exposure_time).
    """

    instrClass = None
    instrDescriptors = ()

    def set_applicable(self):
        self.applicable = []

        if 'PROCESSED_SCIENCE' in self.types:
            return

        if self.descriptors['observation_type'] != 'BPM':
            self.applicable.append('processed_bpm')

        if self.descriptors['observation_type'] != 'DARK':
            self.applicable.append('processed_dark')

        if self.descriptors['observation_type'] not in ('DARK', 'FLAT'):
            self.applicable.append('processed_flat')

    def flat(self, processed=False, howmany=None):
        """Find matching flat: same arm, closest in time."""
        howmany = howmany if howmany else 1

        query = self.get_query() \
            .flat(processed) \
            .match_descriptors(Header.instrument, Header.camera)

        return query.all(howmany)

    def dark(self, processed=False, howmany=None):
        """Find matching dark: same arm, same exposure time, closest in time."""
        howmany = howmany if howmany else 1

        query = self.get_query() \
            .dark(processed) \
            .match_descriptors(Header.instrument, Header.camera,
                               Header.exposure_time)

        return query.all(howmany)

    def bpm(self, processed=False, howmany=None):
        """Find matching BPM: same arm."""
        howmany = howmany if howmany else 1

        query = self.get_query() \
            .bpm(processed) \
            .match_descriptors(Header.instrument, Header.camera)

        return query.all(howmany)

    def dark_coeff(self, processed=False, howmany=None):
        """Find matching dark coefficient file: same arm, closest in time.

        Dark coefficient files are tagged DARK_COEFF.  Like WAVECAL,
        this type has no dedicated FitsStorage REDUCTION_STATUS, so we
        filter by Header.types.
        """
        howmany = howmany if howmany else 1

        query = self.get_query()
        if processed:
            query = query.filter(
                Header.types.contains('PROCESSED'),
                Header.types.contains('DARK_COEFF'),
            )
        else:
            query = query.raw().filter(
                Header.types.contains('DARK_COEFF'),
            )
        query = query.match_descriptors(Header.instrument, Header.camera)

        return query.all(howmany)

    def wavecal(self, processed=False, howmany=None):
        """Find matching wavecal: same arm, closest in time.

        FitsStorage's REDUCTION_STATUS doesn't include WAVECAL, so
        processed wavecals are ingested as PROCESSED_UNKNOWN.  Filter
        by Header.types instead of Header.reduction until the WAVECAL
        tag is migrated to ARC.
        """
        howmany = howmany if howmany else 1

        query = self.get_query()
        if processed:
            query = query.filter(
                Header.types.contains('PROCESSED'),
                Header.types.contains('WAVECAL'),
            )
        else:
            query = query.raw().filter(
                Header.types.contains('WAVECAL'),
            )
        query = query.match_descriptors(Header.instrument, Header.camera)

        return query.all(howmany)
