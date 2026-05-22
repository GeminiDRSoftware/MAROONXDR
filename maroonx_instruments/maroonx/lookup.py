# If the filter names and central wavelength are different from
# the global definitions in gemini_instruments/gemini/lookup.py
# redefine them here in filter_wavelengths.

# all entries here are static
read_noise = {'Q1': 1.14,
    'Q2': 1.14,
    'Q3': 1.14,
    'Q4': 1.14,
    'R1': 1.63,
    'R2': 1.63
} # read noise is actually quoted in variance
gain = {'Q1': 2.72,
    'Q2': 2.72,
    'Q3': 2.72,
    'Q4': 2.72,
    'R1': 2.74,
    'R2': 2.74
}

array_name_b = [
    'Q1', 'Q2', 'Q3', 'Q4'
]

array_name_r = [
    'R1', 'R2'
]
data_section = {  # section of original frame with real pixels
    'Q1': '[7:2040,9:2040]',
    'Q2': '[7:2040,2361:4400]',
    'Q3': '[2361:4280,9:2040]',
    'Q4': '[2361:4280,2361:4400]',
    'R1': '[36:2043,1:4080]',
    'R2': '[2358:4385,1:4080]'
}


array_section = {
    'Q1': '[1:2034,1:2032]',   # where real data will go after overscan removal
    'Q2': '[1:2034,2033:4072]',
    'Q3': '[2035:3954,1:2032]',
    'Q4': '[2035:3954,2033:4072]',
    'R1': '[1:2008,1:4080]',
    'R2': '[2009:4036,1:4080]'
}

# In DRAGONS, bias_section is used for subtracting and trimming the overscan region
# but in the og pipeling, different sections are used for subtracting and trimming
bias_section = {
    'Q1': '[2043:2200,1:2200]',  # overscan location on original frame
    'Q2': '[2043:2200,2201:4400]',
    'Q3': '[2201:2358,1:2200]',
    'Q4': '[2201:2358,2201:4400]',
    'R1': '[2046:2200,1:4400]',
    'R2': '[2201:2355,1:4400]'
}

# Definitions for backward compatibility ============================
# Old section used only for bias subtraction, same for both arms
# array_subtraction_section - mean(bias_subtraction_section)
bias_subtraction_section = {
    'RB1': '[2051:2200,1:2200]',
    'RB2': '[2201:2350,1:2200]',
    'RB3': '[2051:2200,2201:4400]',
    'RB4': '[2201:2350,2201:4400]',
}
array_subtraction_section = {
    'RB1': '[1:2200,1:2200]',
    'RB2': '[2201:4400,1:2200]',
    'RB3': '[1:2200,2201:4400]',
    'RB4': '[2201:4400,2201:4400]',
}
# ===================================================================

detector_section = { 'RB': '[1:4400, 1:4400]'}  #original full-frame, same for both red and blue arms


filter_wavelengths = {
#    'i' : 0.80,
}