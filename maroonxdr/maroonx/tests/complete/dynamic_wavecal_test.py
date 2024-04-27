
import glob
from recipe_system.reduction.coreReduce import Reduce
import sys
import astrodata

from pathlib import Path
parent_dir = Path(__file__).parents[4]
sys.path.append(str(parent_dir))
import maroonx_instruments

all_files = glob.glob('/Users/rohan/Documents/20240411T153824Z_DEEEE_b_0030.fits')

all_files.sort()

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(all_files)
myreduce.drpkg= 'maroonxdr'
myreduce.runr()