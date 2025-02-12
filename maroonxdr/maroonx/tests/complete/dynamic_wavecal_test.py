import glob
import sys
from pathlib import Path

from recipe_system.reduction.coreReduce import Reduce

parent_dir = Path(__file__).parents[4]
sys.path.append(str(parent_dir))


all_files = glob.glob(
    '/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/20220725T031211Z_DEEEE_r_0003_dynamic_wavecal.fits'
)
all_files.sort()

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(all_files)
myreduce.drpkg = 'maroonxdr'
myreduce.runr()
