# If the filter names and central wavelength are different from
# the global definitions in gemini_instruments/gemini/lookup.py
# redefine them here in filter_wavelengths.

# all entries here are static -- !should these stay here or be in headers?
read_noise = [1.1]
gain = [2.7]


array_name_b = [
    'Q1', 'Q2', 'Q3', 'Q4'
]  # put this in headers

array_name_r = [
    'R1', 'R2'
]
data_section = {
    'Q1': '[1:2040,1:2040]',
    'Q2': '[1:2040,2362:4399]',
    'Q3': '[2362:4399,1:2040]',
    'Q4': '[2362:4399,2362:4399]',
    'R1': '[1:2042,1:4080]',
    'R2': '[2359:4399,1:4080]'
}  # put this in headers


array_section = {
    'Q1': '[1:2200,1:2200]',
    'Q2': '[1:2200,2201:4399]',
    'Q3': '[2201:4399,1:2200]',
    'Q4': '[2201:4399,2201:4399]',
    'R1': '[1:2200,1:4399]',
    'R2': '[2201:4399,1:4399]'
}  # put this in headers

bias_section = {
    'Q1': '[1:2200,2051:2200]',
    'Q2': '[1:2200,2201:2350]',
    'Q3': '[2201:4399,2050:2200]',
    'Q4': '[2201:4399,2201:2350]',
    'R1': '[1:2200,4100:4399]',
    'R2': '[2201:4399,4050:4399]'
}  # put this in headers

filter_wavelengths = {
#    'i' : 0.80,
}