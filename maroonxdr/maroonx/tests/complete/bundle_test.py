from pathlib import Path

from gempy.adlibrary import dataselect
from recipe_system.reduction.coreReduce import Reduce

# Get all files in the science_dir.  Change the path here to suit your installation.
science_dir = Path('/home/martin/Documentos/Projects/MAROONXDR/science_dir')
all_files = list(science_dir.glob('*.fits'))
all_files = [str(p) for p in all_files]
all_files.sort()

bundles = dataselect.select_data(all_files, tags=['BUNDLE'])

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(bundles)
myreduce.drpkg = 'maroonxdr'
myreduce.runr()
