from __future__ import division, print_function
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt
import json

import pandas as pd # import data tables from NIST and store with pandas

temperature_pressure_dict {
    (50, 1.2): {
    "density": 5.4,
    "enthalpy": 4.4,
    }
}
isothermal_dict = {
    "50": { #store each temperature for an isotherm.
        "Pressure":     [[1,2,3],'kg'],
        "Density":      [[1,2,3],'kg'],
        "Volume":       [[1,2,3],'kg'],
        "Energy":       [[1,2,3],'kg'],
        "Enthalpy":     [[1,2,3],'kg'],
        "Entropy":      [[1,2,3],'kg'],
        "Gibbs":        [[1,2,3],'kg'],

    },
    "50.1": { #store each temperature for an isotherm.
        "Pressure":     [[1,2,3],'kg'],

    },
}

# --- Command line arguments --- #

with open('data.json', 'w', encoding='utf-8') as file:
     json.dump(isothermal_dict, file, sort_keys=True, ensure_ascii=False, indent=4)

with open('data.json') as file:
    data = json.load(file)
print(data)
print(sys.getsizeof(data)/1000.0, 'kB of Ram')
html = 'https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C74828&Type=IsoTherm&Digits=5&PLow=5&PHigh=100&PInc=.1&T=300&RefState=DEF&TUnit=K&PUnit=bar&DUnit=mol%2Fl&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm'
#
tables = pd.read_html(html)[0]
#
Pressure = (pd.DataFrame(tables[1:],columns=[1])).values.astype(np.float)
#print(np.core.defchararray.strip(Pressure),'[')
print(Pressure.tolist())
# print(type(tables))
print(type(Pressure))


def find_density(p, T, fluid):
    ''' Cache data and only download if we don't have it already '''
    pass

def find_isotherm_internal_energy(pmin, pmax, T, fluid):
    ''' Cache data and only download if we don't have it already '''
    pass

def find_isotherm_entropy(pmin, pmax, T, fluid):
    ''' Cache data and only download if we don't have it already '''
    pass

def find_isotherm_everything(pmin, pmax, T, fluid):
    ''' Cache data and only download if we don't have it already '''
    pass
