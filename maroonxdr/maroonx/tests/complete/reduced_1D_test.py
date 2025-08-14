"""
This script tests the creation of 1D reduced spectra for MAROON-X DEEEE fits files.  
It does not rely on pytest, and does not produce a success or fail output like pytest does.
Instead, if the reduce runs successfully, this will produce a 1D reduced spectra file that can be 
used for dynamic wavelength calibration.  End users should use this test to test their installation.
"""



# from pathlib import Path

# import astrodata

# from recipe_system.reduction.coreReduce import Reduce
# from gempy.adlibrary import dataselect

# import maroonx_instruments  # noqa : important to load adclass tags

# # Get all files in the science_dir.  Change the path here to suit your installation.
# science_dir = Path('/home/martin/Documentos/Projects/MAROONXDR/science_dir')

# all_files = list(science_dir.glob('*_DEEEE_r_*.fits'))
# all_files = [str(f) for f in all_files]


# all_files.sort()
# print(all_files)

# # Run reduce on all selected files
# myreduce = Reduce()
# myreduce.files.extend(all_files)
# myreduce.drpkg= 'maroonxdr'
# myreduce.runr()




import os
from pathlib import Path

from gempy.adlibrary import dataselect
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa : important to load adclass tags

from gempy.utils import logutils
logutils.config(file_name="test_reduction.log", stomp=False)
log = logutils.get_logger(__name__)
log.setLevel("DEBUG")


def test_reduce_wavecal():

    # Get all files in the science_dir.  Change the path here to suit your installation.
    test_path = Path(os.environ.get("MAROONX_DRAGONS_TEST"))
    science_dir = test_path / 'science_dir'

    # Get all files in the science_dir
    all_files = list(Path(science_dir).glob('*.fits'))
    all_files = [str(p) for p in all_files]
    all_files.sort()

    # Change working directory to science_dir
    original_dir = os.getcwd()
    os.chdir(science_dir)
    
    for arm in ['BLUE', 'RED']:
        only_wavecal = dataselect.select_data(all_files, tags=['RAW', 'WAVECAL', arm])

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(only_wavecal)
        myreduce.drpkg = 'maroonxdr'
        myreduce.runr()

    # Restore original working directory
    os.chdir(original_dir)


if __name__ == '__main__':

    test_reduce_wavecal()