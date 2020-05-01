from __future__ import division, print_function
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd # import data tables from NIST and store with pandas

"""
Generate isothermal data from NIST Thermophysical fluid database using Python 3

Copyright (c) 2019 - 2020 Jordan K. Pommerenck

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

usage = '''arguments: dp pmax T FLUID
  pressures are in bar, T in Kelvin

This program downloads data from the (NIST) Thermophysical properties of fluid
system database. This data can then be read using python with pandas or the
users program of choice.

Typical usage:
    python3 isothermal-save-gas-csv.py 0.1 2000 298 xenon
    # The above only needs to be done once
    python3 isothermal-ads-gas-plot.py xenon 298 5.8 200
'''

# --- Command line arguments --- #
if len(sys.argv) != 5:
    print(usage)
    exit(1)
pInc = float(sys.argv[1])
pMax = float(sys.argv[2])
Temperature = float(sys.argv[3])
fluid = sys.argv[4]

basename = '%s-%g' % (fluid, Temperature)
csvname = basename + '.csv'

# The nist_ids dictionary holds the CAS registry identifiers for a few
# interesting species from available species (NIST).
nist_ids = {
    'water':            'C7732185',
    'H20':              'C7732185',
    'nitrogen':         'C7727379',
    'N2':               'C7727379',
    'hydrogen':         'C1333740',
    'H2':               'C1333740',
    'parahydrogen':     'B5000001',
    'deuterium':        '7782390',
    'D2':               '7782390',
    'oxygen':           '7782447',
    'O2':               '7782447',
    'flourine':         'C7782414',
    'F':                'C7782414',
    'carbon-monoxide':  'C630080',
    'CO':               'C630080',
    'carbon-dioxide':   'C124389',
    'CO2':              'C124389',
    'methanol':         'C67561',
    'CH3OH':            'C67561',
    'CH4O':             'C67561',
    'methane':          'C74828',
    'CH4':              'C74828',
    'ethane':           'C74840',
    'C2H6':             'C74840',
    'ethene':           'C74851',
    'propane':          'C74986',
    'C3H8':             'C74986',
    'propene':          'C115071',
    'propyne':          'C74997',
    'cyclopropane':     'C75194',
    'butane':           'C106978',
    'isobutane':        'C75285',
    'pentane':          'C109660',
    'heptane':          'C142825',
    'octane':           'C111659',
    'nonane':           'C111842',
    'decane':           'C124185',
    'dodecane':         'C112403',
    'helium':           'C7440597',
    'neon':             'C7440019',
    'argon':            'C7440371',
    'krypton':          'C7439909',
    'xenon':            'C7440633',
    'ammonia':          'C7664417',
    'benzene':          'C71432',
    'toluene':          'C108883',
    'sulfur-dioxide':   'C7446095',
    'hydrogen-sulfide': 'C7783064',
}

# Check the fluid against the keys in the dictionary. If the input key is not
# found then the user is instructed to add the key to the dictionary.
if (nist_ids.get(fluid)):
    my_id = nist_ids[fluid]
    print(f'{fluid} has a CAS registry of {my_id}.')
else:
    raise KeyError(f"{fluid} is not a valid fluid. You must add name and CAS registry to nist_ids!")

open('data/' + csvname, 'w').close() # empty the file I am creating

num_rows = 500

inc_by = [2, 2, 2, 10/8.0]
inc_counter = 0

p = pInc

while p < pMax:
    print('pressure =', p, 'pressure increment =', pInc)
    p_final = min(pMax, p+pInc*(num_rows-1))
    if p == pInc:
        p_final = min(pMax, p+pInc*(num_rows-2))
    html_web = "https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=%s&" % my_id
    html_arg = "Type=IsoTherm&Digits=12&PLow=%g&PHigh=%g&PInc=%g&T=%g" % (p,
                                                                          p_final,
                                                                          pInc,
                                                                          Temperature)
    html_units = "&RefState=DEF&TUnit=K&PUnit=bar&DUnit=mol%2Fl&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
    # print(html_arg)
    tables = pd.read_html(html_web + html_arg + html_units)[0]
    if p != pInc:
        tables = tables.iloc[1:]
    tables = tables.iloc[:,:7]
    saveData = pd.DataFrame(data = tables)
    saveData.to_csv('data/' + csvname, header=None, index=False, sep=' ', mode='a') # mode = append
    p = p_final + pInc

    pInc *= inc_by[inc_counter % len(inc_by)]
    inc_counter += 1
print('Saved data as', csvname)
