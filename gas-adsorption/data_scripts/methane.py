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

molecular_weight = 16.04*nu.g/mol
density_units = 0.044135*mol/L # density at STP in mol/L
density_unit_name = 'cm$^3$ CH$_4$ STP/cm$^3$'
Name = 'Methane'

n_goals = []
n_arpa_methane_goal = 315*0.044135*mol/L # 315* density at STP
n_goals.append((n_arpa_methane_goal, 'ARPA-E', '#999999', '-'))

def linear_interpolate(x, x1, y1, x2, y2):
    return ((x2-x)*y1 + (x-x1)*y2)/(x2-x1)

# Methane experimental data is from Mason and Long 2014
# The COF-102 data is from https://doi.org/10.1021/ja9015765
# Crys_Dens from DOI: 10.1126/science.1139915
COF102_Crys_Dens = 0.41*gram/cm**3
# I took cm^3/g and converted to L/mol using molecular weight
COF_density_units = mg/gram * COF102_Crys_Dens / molecular_weight

CH4_298_Data = {
    'MOF5': {
        'q_empty': 12.2*kJ/mol,
        'q_full': 12.2*kJ/mol,
        'rho_full': linear_interpolate(65.5, 64.91, 215.14, 69.90, 222.03)*density_units,
        'rho_empty': linear_interpolate(5.8 , 5.58, 28.97, 7.57, 39.00)*density_units
    },
    'HKUST1': {
        'q_empty': 17*kJ/mol,
        'q_full': 17*kJ/mol,
        'rho_full': linear_interpolate(65.5, 65.05, 263.01, 70.01, 266.18)*density_units,
        'rho_empty': linear_interpolate(5.8 , 5.57 , 79.60 ,  7.71, 101.94)*density_units
    },
    'PCN14': {
        'q_empty': 17.6*kJ/mol,
        'q_full': 17.6*kJ/mol,
        'rho_full': linear_interpolate(65.5, 64.89, 239.42, 69.91, 243.67)*density_units,
        'rho_empty': linear_interpolate(5.8 , 5.58, 81.23, 7.56, 100.16)*density_units
    },
    'Ni2dobdc': {
        'q_empty': 20.7*kJ/mol,
        'q_full': 20.7*kJ/mol,
        'rho_full': linear_interpolate(65.5, 65.15, 259.64, 70.06, 260.58)*density_units,
        'rho_empty': linear_interpolate(5.8 , 5.79, 126.15, 7.73 , 144.71)*density_units
    },

    'Co2dobdc': {
        'q_empty': 19.5*kJ/mol,
        'q_full': 19.5*kJ/mol,
        'rho_full': linear_interpolate(65.5, 64.89, 248.64, 69.63, 252.11)*density_units,
        'rho_empty': linear_interpolate(5.8 , 5.62, 117.78, 7.64, 137.04)*density_units
    },
    'Mg2dobdc': {
        'q_empty': 18.5*kJ/mol,
        'q_full': 18.5*kJ/mol,
        'rho_full': linear_interpolate(65.5, 64.90, 230.28, 69.93, 233.71)*density_units,
        'rho_empty': linear_interpolate(5.8 , 5.64, 93.50, 7.60, 111.58)*density_units
    },
    'AX21': {
        'q_empty': 14*kJ/mol,
        'q_full': 14*kJ/mol,
        'rho_full': linear_interpolate(65.5, 64.87, 202.80, 69.65, 209.63)*density_units,
        'rho_empty': linear_interpolate(5.8 , 5.55, 52.19, 7.53, 64.41)*density_units
    },
    # # include the COF-102 data from 10.1021/ja9015765
    # # I used WebPlotDigitizer and extracted the data and interpolated!
    'COF102': {
        'q_empty': 12*kJ/mol,
        'q_full': 12*kJ/mol,
        'rho_full': linear_interpolate(65.5, 59.49, 231.87, 66.03, 237.75)*COF_density_units,
        'rho_empty': linear_interpolate(5.8 , 4.07, 27.65, 6.73, 43.52)*COF_density_units
    },
    # # include the mono-HKUST-1 data from https://doi.org/10.1038/nmat5050
    # # I extracted the data from the supplemental information file TABLE 4.
    'HKUST1-mono': {
        'q_empty': 17*kJ/mol, # same as HKUST-1 ?
        'q_full': 17*kJ/mol, # same as HKUST-1 ?
        'rho_full': linear_interpolate(65.5, 59.79, 254, 69.76, 267)*density_units,
        'rho_empty': linear_interpolate(5.8 , 5.10 , 77 ,  7.10, 98)*density_units
    },
    # # include the MOF-519 data from https://pubs.acs.org/doi/10.1021/ja501606h
    # # I used WebPlotDigitizer on supplemental information file Figure S13.
    'MOF519': {
        'q_empty': 14*kJ/mol, # Extracted from Figure S8
        'q_full': 13*kJ/mol, # Extracted from Figure S8
        'rho_full': linear_interpolate(65.5, 64.65, 261.33, 69.57, 268.28)*density_units,
        'rho_empty': linear_interpolate(5.8 , 3.86 , 39.73 ,  6.91, 65.24)*density_units
    },
}


def isotherm_experiments(T, p_empty_bar, p_full_bar):
    print(T/nu.K)
    if abs(T/nu.K - 298) < 1 and p_empty_bar == 5.8 and p_full_bar == 65:
        return CH4_298_Data
    return {}
