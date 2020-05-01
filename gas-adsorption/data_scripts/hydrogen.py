"""
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
import numpy as np
import matplotlib.pyplot as plt
import numericalunits as nu
from numericalunits import kg, L, mol, mmol, cm, kJ, mg
gram = nu.g

molecular_weight = 2.016*nu.g/nu.mol
density_units = gram/L/molecular_weight
density_unit_name = 'g Hydrogen/L'
Name = 'Hydrogen'

n_goals = []
n_goals.append((0.05*kg/L/molecular_weight, 'DOE ULTIMATE','#999999', '-'))
n_goals.append((0.04*kg/L/molecular_weight, 'DOE 2025',    '#888888', '-'))
n_goals.append((0.03*kg/L/molecular_weight, 'DOE 2020',    '#777777', '-'))

H2_298_Data = {
    'HKUST1': {
        'rho_full': 6.08*mmol/gram,
        'rho_empty': 0.36*mmol/gram,
    },
    'NOTT112': {
        'rho_full': 8.83*mmol/gram,
        'rho_empty': 0.47*mmol/gram,
    },
    'NU125': {
        'rho_full': 8.16*mmol/gram,
        'rho_empty': 0.51*mmol/gram,
    },
    'rhtMOF7': {
        'rho_full': 5.72*mmol/gram,
        'rho_empty': 0.34*mmol/gram,
    },
    'CuMOF74': {
        'rho_full': 3.78*mmol/gram,
        'rho_empty': 0.22*mmol/gram,
    },
    'PCN250': {
        'rho_full': 5.75*mmol/gram,
        'rho_empty': 0.34*mmol/gram,
    },
    'NU1000': {
        'rho_full': 8.67*mmol/gram,
        'rho_empty': 0.51*mmol/gram,
    },
    'UiO67': {
        'rho_full': 5.74*mmol/gram,
        'rho_empty': 0.33*mmol/gram,
    },
    'UiO68Ant': {
        'rho_full': 7.31*mmol/gram,
        'rho_empty': 0.38*mmol/gram,
    },
    'CYCU3Al': {
        'rho_full': 8.91*mmol/gram,
        'rho_empty': 0.49*mmol/gram,
    },
    'Zn2bdc2dabco': {
        'rho_full': 5.29*mmol/gram,
        'rho_empty': 0.31*mmol/gram,
    },
}

H2_77_Data = { # from Table S4. In the supplemental.
    'HKUST1': {
        'q_empty': 4.5*kJ/mol,
        'q_full': 4*kJ/mol,
        'rho_full': 28.04*mmol/gram,
        'rho_empty': 18.32*mmol/gram,
    },
    'NOTT112': {
        'q_empty': 3.6*kJ/mol,
        'q_full': 3.3*kJ/mol,
        'rho_full': 47.86*mmol/gram,
        'rho_empty': 21.51*mmol/gram,
    },
    'NU125': {
        'q_empty': 4.0*kJ/mol,
        'q_full': 3.5*kJ/mol,
        'rho_full': 44.67*mmol/gram,
        'rho_empty': 24.11*mmol/gram,
    },
    'rhtMOF7': {
        'q_empty': 4.4*kJ/mol,
        'q_full': 3.2*kJ/mol,
        'rho_full': 25.68*mmol/gram,
        'rho_empty': 16.71*mmol/gram,
    },
    'CuMOF74': {
        'q_empty': 4.75*kJ/mol,
        'q_full': 4*kJ/mol,
        'rho_full': 16.26*mmol/gram,
        'rho_empty': 11.36*mmol/gram,
    },
    'PCN250': {
        'q_empty': 5.6*kJ/mol,
        'q_full': 4.4*kJ/mol,
        'rho_full': 28.35*mmol/gram,
        'rho_empty': 19.58*mmol/gram,
    },
    'NU1000': {
        'q_empty': 3.4*kJ/mol,
        'q_full': 1.6*kJ/mol,
        'rho_full': 43.38*mmol/gram,
        'rho_empty': 17.46*mmol/gram,
    },
    'UiO67': {
        'q_empty': 4*kJ/mol,
        'q_full': 3*kJ/mol,
        'rho_full': 31.14*mmol/gram,
        'rho_empty': 16.69*mmol/gram,
    },
    'UiO68Ant': {
        'q_empty': 3.75*kJ/mol,
        'q_full': 3*kJ/mol,
        'rho_full': 41*mmol/gram,
        'rho_empty': 19.81*mmol/gram,
    },
    'CYCU3Al': {
        'q_empty': 3.25*kJ/mol,
        'q_full': 2*kJ/mol,
        'rho_full': 44.92*mmol/gram,
        'rho_empty': 17.42*mmol/gram,
    },
    'Zn2bdc2dabco': {
        'q_empty': 4.7*kJ/mol,
        'q_full': 4.5*kJ/mol,
        'rho_full': 25.72*mmol/gram,
        'rho_empty': 17.93*mmol/gram,
    },
    'COF102': {
        'q_empty': 3.5*kJ/mol,  #-> 4.7*kJ/mol,
        'q_full': 3.5*kJ/mol,  # only have low pressure so don't know?
        'rho_full': 62*mg/gram/molecular_weight, #-> 25.72*mmol/gram,
        'rho_empty': 43*mg/gram/molecular_weight, #-> 17.93*mmol/gram,
    },
}

def isotherm_experiments(T, p_empty_bar, p_full_bar):
    print(T/nu.K)
    if abs(T/nu.K - 77) < 1 and p_empty_bar == 5 and p_full_bar == 100:
        return H2_77_Data
    elif abs(T/nu.K - 298) < 1 and p_empty_bar == 5 and p_full_bar == 100:
        return H2_298_Data
    return {}

    
