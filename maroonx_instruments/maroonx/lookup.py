# If the filter names and central wavelength are different from
# the global definitions in gemini_instruments/gemini/lookup.py
# redefine them here in filter_wavelengths.

# data layout and overscan section are static
ccd_section = {
    'Q1': '[0:2200,0:2200]',
    'Q2': '[0:2200,2200:4399]',
    'Q3': '[2200:4399,0:2200]',
    'Q4': '[2200:4399,2200:4399]'
}

bias_section = {
    'Q1': '[0:2200,2050:2200]',
    'Q2': '[0:2200,2200:2350]',
    'Q3': '[2200:4399,2050:2200]',
    'Q4': '[2200:4399,2200:2350]'
}

filter_wavelengths = {
#    'i' : 0.80,
}