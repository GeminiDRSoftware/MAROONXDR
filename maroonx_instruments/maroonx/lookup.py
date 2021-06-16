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
    'Q1': '[1:2049,1:2049]',
    'Q2': '[1:2049,2353:4401]',
    'Q3': '[2353:4401,1:2049]',
    'Q4': '[2353:4401,2353:4401]',
    'R1': '[1:2049,1:4401]',
    'R2': '[2353:4401,1:4401]'
}  # put this in headers


array_section = {
    # 'Q1': '[1:2200,1:2200]',
    # 'Q2': '[1:2200,2201:4399]',
    # 'Q3': '[2201:4399,1:2200]',
    # 'Q4': '[2201:4399,2201:4399]',
    # 'R1': '[1:2200,1:4399]',
    # 'R2': '[2201:4399,1:4399]'
    'Q1': '[1:2049,1:2049]',  # where real data will go after overscan removal
    'Q2': '[1:2049,2049:4097]',
    'Q3': '[2049:4097,1:2049]',
    'Q4': '[2049:4097,2049:4097]',
    'R1': '[1:2049,1:4097]',
    'R2': '[2049:4097,1:4097]'
}  # put this in headers

bias_section = {
    'Q1': '[2051:2201,1:2201]',  # overscan location on original frame
    'Q2': '[2051:2201,2201:2351]',
    'Q3': '[2201:2351,1:2201]',
    'Q4': '[2201:2351,2201:2351]',
    'R1': '[2051:2201,1:4401]',
    'R2': '[2201:2351,1:4401]'
}  # put this in headers

filter_wavelengths = {
#    'i' : 0.80,
}