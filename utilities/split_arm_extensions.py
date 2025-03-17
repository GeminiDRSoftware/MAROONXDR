import copy
from pathlib import Path

import astrodata
from astropy.io import fits


def split_arm_extensions(ad):
    extensions = {}

    for ext in range(len(ad)):
        hdr = ad[ext].hdr
        data = ad[ext].data
        hdu = fits.ImageHDU(data=data, header=hdr, name='SCI')

        # New astrodata for this extension
        arm_ad = astrodata.create(ad.phu, [hdu])

        # Revert to original file name
        arm_name = hdr['ORIGNAME']

        extensions[arm_name] = arm_ad

    return extensions


# This will take all new maroonx files, separate red/blue frames and write
# them back with the original name

base_path = Path('/home/martin/Documentos/Projects/MAROONXDR/science_dir')

subdirs = ['20241114', '20241115']  # , '20241116', '20241124']

for subdir in subdirs:
    files = list((base_path / subdir).glob('*'))

    for file in files:
        ad = astrodata.open(file)

        arm_dict = split_arm_extensions(ad)

        for arm_name, arm_ad in arm_dict.items():
            new_filename = str(base_path / arm_name)

            # Overwrite ORIGNAME from PHU for consistency
            arm_ad.phu['ORIGNAME'] = arm_name

            # Now add a ARCHNAME for reference after ORIGNAME
            archive_card = fits.Card('ARCHNAME', ad.filename, 'Gemini archive name')
            arm_ad.phu.insert('ORIGNAME', archive_card, after=True)

            arm_ad.write(new_filename, overwrite=True)


def splitBundle(adinputs):
    adoutputs = []

    for ad in adinputs:
        for ext in ad.indices:
            hdu = fits.ImageHDU(
                data=ad[ext].data,
                header=ad[ext].hdr,
                name=ad[ext].hdr.get('EXTNAME', 'SCI'),
            )

            # New astrodata for this extension
            arm_ad = astrodata.create(copy.deepcopy(ad.phu), [hdu])

            # Revert to original file name
            arm_name = ad[ext].hdr.get('ORIGNAME')

            # new_filename = str(base_path / arm_name)

            # Overwrite ORIGNAME from PHU for consistency
            arm_ad.phu['ORIGNAME'] = arm_name

            # Now add a ARCHNAME for reference after ORIGNAME
            archive_card = fits.Card('ARCHNAME', ad.filename, 'Gemini archive name')
            arm_ad.phu.insert('ORIGNAME', archive_card, after=True)

            # arm_ad.write(new_filename, overwrite=True)
            adoutputs.append(arm_ad)

    return adoutputs
