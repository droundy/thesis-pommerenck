from __future__ import division, print_function
import sys, os, re, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import numericalunits

import colors

matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')

# --- Build unit system --- #
Kelvin = numericalunits.K
atm = numericalunits.atm
bar = numericalunits.bar
gram = numericalunits.g
kg = numericalunits.kg
L = numericalunits.L
Joule = numericalunits.J
mol = numericalunits.mol
psi = numericalunits.psi
cm = numericalunits.cm
mmol = numericalunits.mmol

mg = numericalunits.mg    # milligram
mL = numericalunits.mL    # milliliter
kJ = numericalunits.kJ     # kiloJoule

# HOW TO RUN:
# python3 gibbs-qst.py

# Define a Global Boltzmann constant
kB = 8.314*Joule/mol/Kelvin

# Define interpolation functions
def my_interp(x0,x,y):
    return np.interp(x0, x, y)

def mu_from_n(particular_n):
    return my_interp(particular_n, n, mu)

def n_from_mu(particular_mu):
    return my_interp(particular_mu, mu, n)

def linear_interpolate(x, x1, y1, x2, y2):
    return ((x2-x)*y1 + (x-x1)*y2)/(x2-x1)

#A list of the five MOFs example AX21-298.txt
#[AX21,Mg2dobdc,Co2dobdc,MOF5,Ni2dobdc,HKUST1,PCN14]

# For the Hydrogen Data (I have to get Isotherms manually using
# online plot digitizer although I can save the online project file
# as data)

#rho_MOF = 18.32*mmol/gram*molecular_weight/MOF_DENSITY

# Looks like I need Qst at p_empty and p_full since most of the curves
# are not constant!
H2_crystal_density = {
    'HKUST1': 0.878*gram/cm**3,
    'NOTT112': 0.446*gram/cm**3,
    'NU125': 0.578*gram/cm**3,
    'rhtMOF7': 0.789*gram/cm**3,
    'CuMOF74': 01.323*gram/cm**3,
    'PCN250': 0.896*gram/cm**3,
    'NU1000': 0.571*gram/cm**3,
    'UiO67': 0.688*gram/cm**3,
    'UiO68Ant': 0.607*gram/cm**3,
    'CYCU3Al': 0.447*gram/cm**3,
    'Zn2bdc2dabco': 0.873*gram/cm**3,
    'NU1101': 0.459*gram/cm**3,
    'NU1102': 0.403*gram/cm**3,
    'NU1103': 0.298*gram/cm**3,
}

# H2_77_160_Data = {
    # 'HKUST1': {
        # 'qh': 6.3*kJ/mol,
        # 'qc': 4*kJ/mol,
        # 'rho_c': 28.04*mmol/gram,
        # 'rho_h': 2.23*mmol/gram,
    # },
    # 'NOTT112': {
        # 'qh': 4.8*kJ/mol,
        # 'qc': 3.3*kJ/mol,
        # 'rho_c': 47.86*mmol/gram,
        # 'rho_h': 2.54*mmol/gram,
    # },
    # 'NU125': {
        # 'qh': 4.8*kJ/mol,
        # 'qc': 3.5*kJ/mol,
        # 'rho_c': 44.67*mmol/gram,
        # 'rho_h': 2.62*mmol/gram,
    # },
    # 'rhtMOF7': {
        # 'qh': 5.9*kJ/mol,
        # 'qc': 3.2*kJ/mol,
        # 'rho_c': 25.68*mmol/gram,
        # 'rho_h': 2.13*mmol/gram,
    # },
    # 'CuMOF74': {
        # 'qh': 5.8*kJ/mol,
        # 'qc': 4*kJ/mol,
        # 'rho_c': 16.26*mmol/gram,
        # 'rho_h': 1.51*mmol/gram,
    # },
    # 'PCN250': {
        # 'qh': 6.2*kJ/mol,
        # 'qc': 4.4*kJ/mol,
        # 'rho_c': 28.35*mmol/gram,
        # 'rho_h': 2.43*mmol/gram,
    # },
    # 'NU1000': {
        # 'qh': 4.7*kJ/mol,
        # 'qc': 1.6*kJ/mol,
        # 'rho_c': 43.38*mmol/gram,
        # 'rho_h': 2.115*mmol/gram,
    # },
    # 'UiO67': {
        # 'qh': 5.5*kJ/mol,
        # 'qc': 3*kJ/mol,
        # 'rho_c': 31.14*mmol/gram,
        # 'rho_h': 1.85*mmol/gram,
    # },
    # 'UiO68Ant': {
        # 'qh': 5.5*kJ/mol,
        # 'qc': 3*kJ/mol,
        # 'rho_c': 41*mmol/gram,
        # 'rho_h': 2.34*mmol/gram,
    # },
    # 'CYCU3Al': {
        # 'qh': 4.4*kJ/mol,
        # 'qc': 2*kJ/mol,
        # 'rho_c': 44.92*mmol/gram,
        # 'rho_h': 1.90*mmol/gram,
    # },
    # 'Zn2bdc2dabco': {
        # 'qh': 4.9*kJ/mol,
        # 'qc': 4.5*kJ/mol,
        # 'rho_c': 25.72*mmol/gram,
        # 'rho_h': 1.65*mmol/gram,
    # },
# }

# from Table S4. In the supplemental. 
# https://doi.org/10.1021/acsenergylett.8b00154

H2_298_Data = {
    'HKUST1': {
        'qh': 4.5*kJ/mol,
        'qc': 4*kJ/mol,
        'rho_c': 6.08*mmol/gram,
        'rho_h': 0.36*mmol/gram,
    },
    'NOTT112': {
        'qh': 3.6*kJ/mol,
        'qc': 3.3*kJ/mol,
        'rho_c': 8.83*mmol/gram,
        'rho_h': 0.47*mmol/gram,
    },
    'NU125': {
        'qh': 4.0*kJ/mol,
        'qc': 3.5*kJ/mol,
        'rho_c': 8.16*mmol/gram,
        'rho_h': 0.51*mmol/gram,
    },
    'rhtMOF7': {
        'qh': 4.4*kJ/mol,
        'qc': 3.2*kJ/mol,
        'rho_c': 5.72*mmol/gram,
        'rho_h': 0.34*mmol/gram,
    },
    'CuMOF74': {
        'qh': 4.75*kJ/mol,
        'qc': 4*kJ/mol,
        'rho_c': 3.78*mmol/gram,
        'rho_h': 0.22*mmol/gram,
    },
    'PCN250': {
        'qh': 5.6*kJ/mol,
        'qc': 4.4*kJ/mol,
        'rho_c': 5.75*mmol/gram,
        'rho_h': 0.34*mmol/gram,
    },
    'NU1000': {
        'qh': 3.4*kJ/mol,
        'qc': 1.6*kJ/mol,
        'rho_c': 8.67*mmol/gram,
        'rho_h': 0.51*mmol/gram,
    },
    'UiO67': {
        'qh': 4*kJ/mol,
        'qc': 3*kJ/mol,
        'rho_c': 5.74*mmol/gram,
        'rho_h': 0.33*mmol/gram,
    },
    'UiO68Ant': {
        'qh': 3.75*kJ/mol,
        'qc': 3*kJ/mol,
        'rho_c': 7.31*mmol/gram,
        'rho_h': 0.38*mmol/gram,
    },
    'CYCU3Al': {
        'qh': 3.25*kJ/mol,
        'qc': 2*kJ/mol,
        'rho_c': 8.91*mmol/gram,
        'rho_h': 0.49*mmol/gram,
    },
    'Zn2bdc2dabco': {
        'qh': 4.7*kJ/mol,
        'qc': 4.5*kJ/mol,
        'rho_c': 5.29*mmol/gram,
        'rho_h': 0.31*mmol/gram,
    },
}

# from Table S4. In the supplemental. 
# https://doi.org/10.1021/acsenergylett.8b00154

H2_77_Data = {
    'HKUST1': {
        'qh': 4.5*kJ/mol,
        'qc': 4*kJ/mol,
        'rho_c': 28.04*mmol/gram,
        'rho_h': 18.32*mmol/gram,
    },
    'NOTT112': {
        'qh': 3.6*kJ/mol,
        'qc': 3.3*kJ/mol,
        'rho_c': 47.86*mmol/gram,
        'rho_h': 21.51*mmol/gram,
    },
    'NU125': {
        'qh': 4.0*kJ/mol,
        'qc': 3.5*kJ/mol,
        'rho_c': 44.67*mmol/gram,
        'rho_h': 24.11*mmol/gram,
    },
    'rhtMOF7': {
        'qh': 4.4*kJ/mol,
        'qc': 3.2*kJ/mol,
        'rho_c': 25.68*mmol/gram,
        'rho_h': 16.71*mmol/gram,
    },
    'CuMOF74': {
        'qh': 4.75*kJ/mol,
        'qc': 4*kJ/mol,
        'rho_c': 16.26*mmol/gram,
        'rho_h': 11.36*mmol/gram,
    },
    'PCN250': {
        'qh': 5.6*kJ/mol,
        'qc': 4.4*kJ/mol,
        'rho_c': 28.35*mmol/gram,
        'rho_h': 19.58*mmol/gram,
    },
    'NU1000': {
        'qh': 3.4*kJ/mol,
        'qc': 1.6*kJ/mol,
        'rho_c': 43.38*mmol/gram,
        'rho_h': 17.46*mmol/gram,
    },
    'UiO67': {
        'qh': 4*kJ/mol,
        'qc': 3*kJ/mol,
        'rho_c': 31.14*mmol/gram,
        'rho_h': 16.69*mmol/gram,
    },
    'UiO68Ant': {
        'qh': 3.75*kJ/mol,
        'qc': 3*kJ/mol,
        'rho_c': 41*mmol/gram,
        'rho_h': 19.81*mmol/gram,
    },
    'CYCU3Al': {
        'qh': 3.25*kJ/mol,
        'qc': 2*kJ/mol,
        'rho_c': 44.92*mmol/gram,
        'rho_h': 17.42*mmol/gram,
    },
    'Zn2bdc2dabco': {
        'qh': 4.7*kJ/mol,
        'qc': 4.5*kJ/mol,
        'rho_c': 25.72*mmol/gram,
        'rho_h': 17.93*mmol/gram,
    },
}

# Methane Experimental Data from Mason and Long 2014
CH4_298_E_Data = {
    'HKUST1': 17*kJ/mol,
    'PCN14': 17.6*kJ/mol,
    'Ni2dobdc': 20.7*kJ/mol,
    'Co2dobdc': 19.5*kJ/mol,
    'Mg2dobdc':18.5*kJ/mol,
    'MOF5': 12.2*kJ/mol,
    'AX21': 14*kJ/mol,
}

# --- Read in the data from NIST --- #
NistData = np.loadtxt('data/methane-298.csv', skiprows=1)
molecular_weight = 16.04*gram/mol
density_units = 0.044135*mol/L # density at STP in mol/L
p_empty = 5.8*bar
p_full = 65*bar

T = NistData[:,0]*Kelvin            # Temperature
p = NistData[:,1]*bar               # Pressure
V = NistData[:,3]*L/mol             # Volume
H = NistData[:,5]*kJ/mol            # Enthalpy
S = NistData[:,6]*Joule/mol/Kelvin  # Entropy

# Calculate the free energies
G = H - T*S

# Calculate chemical potential, number, inverse beta
mu = G # per mol
n = 1/V
mu_empty = my_interp(p_empty, p, mu)
mu_full = my_interp(p_full, p, mu)

# --- Read in the data from MOF file --- #

# Note that all of the *-298.txt files are for methane
# total adsorption and come from the Mason and Long 2014 paper
mofs_shown = set()
plt.figure(figsize=(5,5))
for mof in colors.order(CH4_298_E_Data): # For each MOF
    fname = 'data/%s/%s-%g.txt' % ('methane', mof, 298)
    MofData = np.loadtxt(fname, skiprows=1,unpack=True)
    T_MOF = MofData[0][0]*Kelvin
    p_MOF = MofData[1]*bar
    rho_MOF = MofData[2]*density_units
    rho_lo_p = my_interp(p_empty, p_MOF, rho_MOF)
    delta_G_lo_p = mu_from_n(rho_lo_p) - mu_empty

    rho_hi_p = my_interp(p_full, p_MOF, rho_MOF)
    delta_G_hi_p = mu_from_n(rho_hi_p) - mu_full

    qst = CH4_298_E_Data[mof]
    plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
             [qst/(kJ/mol), qst/(kJ/mol)],
             '.-', label=colors.latex_mof(mof), color=colors.color(mof))
    mofs_shown.add(colors.latex_mof(mof))




# --- Read in the data from NIST --- #
NistData = None
mu = None
p = None
NistDataCold = np.loadtxt('data/' + 'hydrogen-77.csv', skiprows=1)
NistDataHot = np.loadtxt('data/' + 'hydrogen-77.csv', skiprows=1)
molecular_weight = 2.016*gram/mol
density_units = kg/L/molecular_weight
p_empty = 5*bar
p_full = 100*bar

Tc = NistDataCold[:,0]*Kelvin            # Temperature
pc = NistDataCold[:,1]*bar               # Pressure
Vc = NistDataCold[:,3]*L/mol             # Volume
Hc = NistDataCold[:,5]*kJ/mol            # Enthalpy
Sc = NistDataCold[:,6]*Joule/mol/Kelvin  # Entropy
Gc = Hc - Tc*Sc
mu_c = Gc # per mol
nc = 1/Vc

Th = NistDataHot[:,0]*Kelvin            # Temperature
ph = NistDataHot[:,1]*bar               # Pressure
Vh = NistDataHot[:,3]*L/mol             # Volume
Hh = NistDataHot[:,5]*kJ/mol            # Enthalpy
Sh = NistDataHot[:,6]*Joule/mol/Kelvin  # Entropy
Gh = Hh - Th*Sh
mu_h = Gh # per mol
nh = 1/Vh

mu_empty = my_interp(p_empty, ph, mu_h)
mu_full = my_interp(p_full, pc, mu_c)

for mof in colors.order(H2_77_Data): # For each MOF
    rho_lo_p = H2_77_Data[mof]['rho_h']*H2_crystal_density[mof]
    delta_G_lo_p = my_interp(rho_lo_p, nh, mu_h) - mu_empty
    # convert units: mmol/gram*molecular_weight/crystal_density,
    rho_hi_p = H2_77_Data[mof]['rho_c']*H2_crystal_density[mof]
    delta_G_hi_p = my_interp(rho_hi_p, nc, mu_c) - mu_full

    qst = H2_77_Data[mof]
    if colors.latex_mof(mof) not in mofs_shown:
        plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
                 [qst['qh']/(kJ/mol), qst['qc']/(kJ/mol)],
                 'x-', label=colors.latex_mof(mof), color=colors.color(mof))
        mofs_shown.add(colors.latex_mof(mof))
    else:
        plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
                 [qst['qh']/(kJ/mol), qst['qc']/(kJ/mol)],
                 'x-', color=colors.color(mof))

NistData = None
mu = None
p = None
NistDataCold = np.loadtxt('data/' + 'hydrogen-298.csv', skiprows=1)
NistDataHot = np.loadtxt('data/' + 'hydrogen-298.csv', skiprows=1)
molecular_weight = 2.016*gram/mol
density_units = kg/L/molecular_weight
p_empty = 5*bar
p_full = 100*bar

Tc = NistDataCold[:,0]*Kelvin            # Temperature
pc = NistDataCold[:,1]*bar               # Pressure
Vc = NistDataCold[:,3]*L/mol             # Volume
Hc = NistDataCold[:,5]*kJ/mol            # Enthalpy
Sc = NistDataCold[:,6]*Joule/mol/Kelvin  # Entropy
Gc = Hc - Tc*Sc
mu_c = Gc # per mol
nc = 1/Vc

Th = NistDataHot[:,0]*Kelvin            # Temperature
ph = NistDataHot[:,1]*bar               # Pressure
Vh = NistDataHot[:,3]*L/mol             # Volume
Hh = NistDataHot[:,5]*kJ/mol            # Enthalpy
Sh = NistDataHot[:,6]*Joule/mol/Kelvin  # Entropy
Gh = Hh - Th*Sh
mu_h = Gh # per mol
nh = 1/Vh

mu_empty = my_interp(p_empty, ph, mu_h)
mu_full = my_interp(p_full, pc, mu_c)

for mof in colors.order(H2_298_Data): # For each MOF
    rho_lo_p = H2_298_Data[mof]['rho_h']*H2_crystal_density[mof]
    delta_G_lo_p = my_interp(rho_lo_p, nh, mu_h) - mu_empty
    # convert units: mmol/gram*molecular_weight/crystal_density,
    rho_hi_p = H2_298_Data[mof]['rho_c']*H2_crystal_density[mof]
    delta_G_hi_p = my_interp(rho_hi_p, nc, mu_c) - mu_full

    qst = H2_298_Data[mof]
    if colors.latex_mof(mof) not in mofs_shown:
        plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
                 [qst['qh']/(kJ/mol), qst['qc']/(kJ/mol)],
                 'x-', label=colors.latex_mof(mof), color=colors.color(mof))
    else:
        plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
                 [qst['qh']/(kJ/mol), qst['qc']/(kJ/mol)],
                 'x-', color=colors.color(mof))

plt.plot([0,31],[0,31],'k:')
plt.xlim(0,22)
plt.ylim(0,22)
plt.axes().set_aspect('equal')
plt.xlabel(r'$|\Delta G_{st}|$ (kJ/mol)')
plt.ylabel(r'$|q_{st}|$ (kJ/mol)')
plt.legend(loc='best', prop=colors.small_font)
plt.savefig('figs/qst-vs-delta-G.pdf')

# matplotlib argument for showing plots
if 'noshow' not in sys.argv:
    plt.show()
