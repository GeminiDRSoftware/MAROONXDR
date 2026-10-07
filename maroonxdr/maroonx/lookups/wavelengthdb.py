# maroonx/wavelengthdb.py
#
# This file contains the reference and static wavelength lookup tables
# for MaroonX. Keys are the arm name ('BLUE' or 'RED'); values are file
# names relative to lookups/WLS/, or absolute paths.

refwavelength_dict = {
    "BLUE" : "REFWAVELENGTH_b.fits",
    "RED" : "REFWAVELENGTH_r.fits"
}

statwavelength_dict = {
    "BLUE" : "WLSTAT_b.fits",
    "RED" : "WLSTAT_r.fits"
}
