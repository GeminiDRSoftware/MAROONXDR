
import glob
from recipe_system.reduction.coreReduce import Reduce
import sys
import astrodata

from pathlib import Path
parent_dir = Path(__file__).parents[4]
sys.path.append(str(parent_dir))
import maroonx_instruments

all_files = glob.glob('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/science_dir/20220725T031211Z_DEEEE_r_0003.fits')

all_files.sort()

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(all_files)
myreduce.drpkg= 'maroonxdr'
myreduce.runr()