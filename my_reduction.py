
from pathlib import Path
import itertools as it

from gempy.adlibrary import dataselect
from recipe_system.reduction.coreReduce import Reduce

import astrodata
import maroonx_instruments  # noqa : important to load adclass tags

from gempy.utils import logutils
logutils.config(file_name="test.log", mode="debug", stomp=True)

# Get all files in the science_dir.  Change the path here to suit your installation.
science_dir = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir')
def get_files(path=None):
    if path is None:
        path = science_dir
    all_files = list(path.glob('*.fits'))
    all_files = [str(f) for f in all_files]
    all_files.sort()
    return all_files

# some testing files from previous runs are in science_dir
# so we need to filter them out
raw_files = dataselect.select_data(get_files(), tags=['RAW'])



# =============================================================================
# Step 1 - Debundle the data
# Select bundles
selected_bundles = dataselect.select_data(get_files(), tags=['BUNDLE'])

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(selected_bundles)
myreduce.drpkg = 'maroonxdr'
myreduce.runr()

# Step 2 - Master Flats
for arm in ['RED', 'BLUE']:

    # Select both DFFFD and FDDDF files
    selected_flats = dataselect.select_data(get_files(), tags=['RAW', 'FLAT', arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_flats)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'makeProcessedFlatDFFFF'
    myreduce.runr()


# Step 3 - Master Darks
exptime_tags = ['300s']
#exptime_tags = ['60s', '120s', '300s', '600s', '900s', '1200s', '1800s']
arm_tags = ['BLUE', 'RED']

for exptime, arm in it.product(exptime_tags, arm_tags):

    selected_darks = dataselect.select_data(get_files(), tags=['DARK', exptime, arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_darks)
    myreduce.drpkg = 'maroonxdr'
    myreduce.runr()


# Step 4 - Extract flux
arm_tags = ['RED', 'BLUE']

for arm in arm_tags:

    selected_spect = dataselect.select_data(get_files(), tags=['WAVECAL', arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_spect)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'makeDynamicWavecal'
    myreduce.runr()


# =============================================================================
# FLAT reduction
# =============================================================================
# Reduction fragments for testing on ipython
import astrodata
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX
from copy import deepcopy

selected_spect = dataselect.select_data(get_files(), tags=['RAW', 'FLAT', 'BLUE'])

# read files and instantiate the primitive class
adinput = [astrodata.open(f) for f in selected_spect]
p = MAROONX(adinput)

p.prepare()
p.checkArm()
p.checkND()
p.addDQ()
p.subtractOverscan()
p.trimOverscan()
p.correctImageOrientation()
p.addVAR(read_noise=True,poisson_noise=True)
p.separateFlatStreams()  # creates 'DFFFD_flats' stream and leaves FDDDF flats in main stream

# stack each stream separately
p.stackFlats(stream='main', scale_mode='mean_frame', suffix='_DDDDF_flat')
p.stackFlats(stream='DFFFD_flats', scale_mode='mean_frame', suffix='_DFFFD_flat')


p.findStripes()  # define stripe info to ultimately remove stray light in each stream
p.findStripes(stream='DFFFD_flats')

p.identifyStripes(selected_fibers='0,0,0,0,5') # identify stripes based on MX architecture files
p.identifyStripes(stream='DFFFD_flats',selected_fibers='0,2,3,4,0')

p.defineFlatStripes()  # defines pixel inclusion for each flat region based on stripe ids
p.defineFlatStripes(stream='DFFFD_flats')

p.removeStrayLight(filter_size=21, box_size=21, snapshot=True)  # remove straylight from frame (this is why 2 partial illumination flat sets are necessary)
p.removeStrayLight(stream='DFFFD_flats', filter_size=21, box_size=21, snapshot=True)
# p.writeOutputs(suffix='_remove_straylight')  # save reduced flat frames

# p.combineFlatStreams(stream='main', stream_2='DFFFD_flats')  # combine straylight-removed images
# p.clearStream(stream='DFFFD_flats') # remove second stream

# p.findStripes()  # re-run find/identify/define routine on combined frame
# p.identifyStripes(selected_fibers='0,2,3,4,5')
# p.defineFlatStripes(extract=True)
# # run the 5-illuminated-fiber frame through extraction to create a reduced processed flat
# p.storeProcessedFlat(suffix='_DFFFF_flat')


# ad = deepcopy(p.streams['main'][0])


# =============================================================================
# Spectrum reduction
# =============================================================================
# Reduction fragments for testing on ipython

from maroonxdr.maroonx import parameters_maroonx_spectrum
from maroonxdr.maroonx.maroonx_fit import maroonx_fit
from maroonxdr.maroonx import maroonx_utils
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle
from maroonxdr.maroonx.maroonx_echellespectrum.maroonxspectrum import MXSpectrum
from maroonxdr.maroonx.maroonx_echellespectrum.wavelengthsolution import WavelengthSolution

from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from copy import deepcopy


selected_spect = dataselect.select_data(get_files(), tags=['RAW', 'WAVECAL', 'BLUE'])
#selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'WAVECAL', 'RED'], xtags=['LFC'])

# selected_spect = dataselect.select_data(get_files(), tags=['PROCESSED'])  # all_files used

# read files and instantiate the primitive class
adinput = [astrodata.open(f) for f in selected_spect[3:6]]
p = MaroonXSpectrum(adinput)

p.prepare()
p.checkArm()
p.addDQ()  # just placeholder until MX is in caldb
p.overscanCorrect()
p.correctImageOrientation()
p.addVAR(read_noise=True, poisson_noise=True)
# # get and save wavelength solution (either static reference or frame's unique sim cal solved)
# first perform echelle extraction of fibers
p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra
p.boxExtraction() # extracts spectra from stripes
p.getPeaksAndPolynomials(fibers=(2,), multithreading=False) # fits etalon peaks and polynomials
# p.writeOutputs(suffix='_dynamic_wavecal')  # save reduced 1D spectra

p.staticWavelengthSolution()
p.fitAndApplyEtalonWls()


# =============================================================================
# Spectrum reduction
# =============================================================================
# Reduction fragments for testing on ipython

from maroonxdr.maroonx import parameters_maroonx_spectrum
from maroonxdr.maroonx.maroonx_fit import maroonx_fit
from maroonxdr.maroonx import maroonx_utils
from maroonxdr.maroonx.primitives_maroonx_echelle import MAROONXEchelle
from maroonxdr.maroonx.maroonx_echellespectrum.maroonxspectrum import MXSpectrum
from maroonxdr.maroonx.maroonx_echellespectrum.wavelengthsolution import WavelengthSolution

from maroonxdr.maroonx.primitives_maroonx_spectrum import MaroonXSpectrum
from copy import deepcopy


selected_spect = dataselect.select_data(get_files(), tags=['RAW', 'SCI', 'BLUE', '300s'])
#selected_spect = dataselect.select_data(raw_files, tags=['RAW', 'WAVECAL', 'RED'], xtags=['LFC'])

# selected_spect = dataselect.select_data(get_files(), tags=['PROCESSED'])  # all_files used

# read files and instantiate the primitive class
adinput = [astrodata.open(f) for f in selected_spect]
p = MaroonXSpectrum(adinput)

p.prepare()
p.checkArm()
p.addDQ()  # just placeholder until MX is in caldb
p.overscanCorrect()
p.correctImageOrientation()
p.addVAR(read_noise=True,poisson_noise=True)
# get and save wavelength solution (either static reference or frame's unique sim cal solved)
p.darkSubtraction()
p.extractStripes()  # gets relevant flat and dark to cut out frame's spectra TODO Skip dark for fiber 5
p.optimalExtraction()  # does 2D to 1D conversion of cut out spectra (only for fibers 2,3,4)
# TODO: perform echelle peak fitting on fiber 5
# TODO: Get wavelength solution from dynamic wavecal recipe
# TODO: Take Fiber 5 peak positions and 
p.storeProcessedScience(suffix='_reduced')