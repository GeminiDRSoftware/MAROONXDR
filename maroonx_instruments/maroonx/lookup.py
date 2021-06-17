# If the filter names and central wavelength are different from
# the global definitions in gemini_instruments/gemini/lookup.py
# redefine them here in filter_wavelengths.

# all entries here are static -- !should these stay here or be in headers?
read_noise = {'Q1': 1.1,
    'Q2': 1.1,
    'Q3': 1.1,
    'Q4': 1.1,
    'R1': 1.1,
    'R2': 1.1
}
gain = {'Q1': 2.7,
    'Q2': 2.7,
    'Q3': 2.7,
    'Q4': 2.7,
    'R1': 2.7,
    'R2': 2.7
}  # not correct as written just hacking for now

array_name_b = [
    'Q1', 'Q2', 'Q3', 'Q4'
]  # put this in headers

array_name_r = [
    'R1', 'R2'
]
data_section = {  # section of original frame with real pixels
    'Q1': '[1:2040,1:2040]',
    'Q2': '[1:2040,2361:4400]',
    'Q3': '[2361:4400,1:2040]',
    'Q4': '[2361:4400,2361:4400]',
    'R1': '[4:2043,1:4080]',
    'R2': '[2358:4397,1:4080]'
}  # put this in headers


array_section = {
    'Q1': '[1:2048,1:2048]',  # where real data will go after overscan removal
    'Q2': '[1:2048,2049:4096]',
    'Q3': '[2049:4096,1:2048]',
    'Q4': '[2049:4096,2049:4096]',
    'R1': '[1:2048,1:4096]',
    'R2': '[2049:4096,1:4096]'
}  # put this in headers

bias_section = {
    'Q1': '[2043:2200,1:2200]',  # overscan location on original frame
    'Q2': '[2043:2200,2201:4400]',
    'Q3': '[2201:2358,1:2200]',
    'Q4': '[2201:2358,2201:4400]',
    'R1': '[2046:2200,1:4400]',
    'R2': '[2201:2355,1:4400]'
}  # put this in headers

filter_wavelengths = {
#    'i' : 0.80,
}