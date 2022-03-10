from __future__ import division, print_function
import sys, os, re, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import numericalunits, colors

matplotlib.rcParams['text.usetex'] = True

from matplotlib.font_manager import FontProperties

small_font = FontProperties()
small_font.set_size('small')

"""
Create deliverable capacity and density plots from NIST Thermophysical fluid
data generated with isothermal-save-gas-csv.py using Python 3

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

# --- Build unit system --- #

# We use the numericalunits package to ensure that our units are consistent.
# This package works by creating randomized units for internal use, such that
# any inconsistency in use of units will produce different output on each
# invocation.
Kelvin = numericalunits.K
atm = numericalunits.atm
bar = numericalunits.bar
gram = numericalunits.g
kg = numericalunits.kg
L = numericalunits.L
Joule = numericalunits.J
mol = 6.02214e23
cm = numericalunits.cm
angstrom = 1e-8*cm
mmol = 1e-3*mol

mg = numericalunits.mg    # milligram
mL = numericalunits.mL    # milliliter
kJ = numericalunits.kJ     # kiloJoule
J = numericalunits.J

# --- Command line arguments --- #

basename = 'methane'
temperature = '298'
p_empty = 5.8*bar
p_full = 65*bar

kB = 8.314*Joule/mol/Kelvin

gas_aliases = {
    'H2O':              'water',
    'N2':               'nitrogen',
    'H2':               'hydrogen',
    'D2':               'deuterium',
    'O2':               'oxygen',
}
if basename in gas_aliases:
    basename = gas_aliases[basename]

# Add the data_scripts directory to the path to import
sys.path.insert(0,'data_scripts/')
gas = __import__(basename)

molecular_weight = gas.molecular_weight

# Choose the appropriate density units.  For methane it is common to plot output
# in volume/volume STP.  For other systems, we just output density in g/L.
density_units = gas.density_units
density_unit_name = gas.density_unit_name

# Load the data from files created by isothermal-save-gas-csv.py
data = np.loadtxt('data/' + '%s-%s.csv' % (basename, temperature), skiprows=1)

# Create thermodynamic variables for pure gas in our fancy units.
T = data[:,0]*Kelvin            # Temperature
p = data[:,1]*bar               # Pressure
rho = data[:,2]*mol/L           # Density
V = data[:,3]*L/mol             # Volume
U = data[:,4]*kJ/mol            # Internal Energy
H = data[:,5]*kJ/mol            # Enthalpy
S = data[:,6]*Joule/mol/Kelvin  # Entropy

G = H - T*S
F = U - T*S

mu = G # per mol
n = 1/V
kT = T*kB

# Define interpolation functions
def my_interp(x, xp, fp):
    interp = np.interp(x, xp, fp, left=np.nan, right=np.nan)
    interp = interp[~np.isnan(interp)]
    if np.isnan(interp).any():
        raise Exception('Trying to evaluate outside the data will crash interpolation!')
    return interp

def mu_from_n(particular_n):
    return my_interp(particular_n, n, mu)

def n_from_mu(particular_mu):
    return my_interp(particular_mu, mu, n)

def mu_from_p(particular_p):
    return my_interp(particular_p, p, mu)

mu_empty = mu_from_p(p_empty)
mu_full = mu_from_p(p_full)

# Initialize max_Gads for setting plot limit.
max_Gads = 1000
Gads = np.linspace(0,max_Gads,10000)*kJ/mol

#------------------------------------------------------------------------------#
# Create the first figure of deliverable capacity vs Detla Gst.
plt.figure('deliverable capacity', figsize=(6,3))
rho_f = n_from_mu(mu_full+Gads)
rho_e = n_from_mu(mu_empty+Gads)[0:len(rho_f)]
D = rho_f - rho_e

# Determine the X and Y limits for the plot
y = D/density_units
ymax = max(y)
xpos = np.where(y == ymax)[0]
x = Gads/(kJ/mol)
xmax = x[xpos][0]

XX = Gads[0:len(D)]/(kJ/mol)
YY = D/density_units

max_index = np.argwhere(YY == np.amax(YY))[0][0]
derivative = np.diff(YY)[max_index:]/np.diff(XX)[max_index:]
minval = np.nanmin(derivative[np.isfinite(derivative)])
new_index = np.argwhere(derivative != 0)[-10][0]

biggest_x_reasonable = Gads[0:new_index].max()/(kJ/mol)
x_max_lim = min(3*xmax, biggest_x_reasonable)
plt.xlim(0,x_max_lim)
y_max_lim = np.max(n_from_mu(mu_full+Gads))/density_units
y_max_lim = 1.5*ymax

bbox = dict(boxstyle="round", fc="0.8")

x_goal_label = 0.03*x_max_lim
for n_goal, label, color, style in gas.n_goals:
    plt.text(x_goal_label, n_goal/density_units+0.02*ymax, label+' target', color='k',
             bbox=dict(facecolor='white', edgecolor='white', alpha=1, pad=0))
    line = plt.axhline(n_goal/density_units, color='k', linestyle=style,
                       linewidth=0.5)

# Create the legend for the plot
to_be_legended = []
legend_labels = []

plt.plot(XX, D/density_units)

# DATA is from TABLE 3 DOI: 10.1021/acsenergylett.8b00154
# Crys_Dens for COF from DOI: 10.1126/science.1139915
crystal_density = {
    'HKUST1':           0.878*gram/cm**3,
    'NOTT112':          0.446*gram/cm**3,
    'NU125':            0.578*gram/cm**3,
    'rhtMOF7':          0.789*gram/cm**3,
    'CuMOF74':          1.323*gram/cm**3,
    'PCN250':           0.896*gram/cm**3,
    'NU1000':           0.571*gram/cm**3,
    'UiO67':            0.688*gram/cm**3,
    'UiO68Ant':         0.607*gram/cm**3,
    'CYCU3Al':          0.447*gram/cm**3,
    'Zn2bdc2dabco':     0.873*gram/cm**3,
    'NU1101':           0.459*gram/cm**3,
    'NU1102':           0.403*gram/cm**3,
    'NU1103':           0.298*gram/cm**3,
    'COF102':           0.41*gram/cm**3,
    # WARNING, the following is false, and is intended to compensate for rho
    # being given in mmol/L not mmol/g
    'Ni2m-dobdc':       1,
}

mof_isotherms = gas.isotherm_experiments(T[0], 5.8, 65)
for mof in colors.order(mof_isotherms): # For each MOF
    if basename == 'methane':
        rho_lo_p = mof_isotherms[mof]['rho_empty']#*crystal_density[mof]
    else:
        rho_lo_p = mof_isotherms[mof]['rho_empty']*crystal_density[mof]
    delta_G_lo_p = np.interp(rho_lo_p, n, mu) - mu_empty

    if basename == 'methane':
        rho_hi_p = mof_isotherms[mof]['rho_full']#*crystal_density[mof]
    else:
        rho_hi_p = mof_isotherms[mof]['rho_full']*crystal_density[mof]
    delta_G_hi_p = np.interp(rho_hi_p, n, mu) - mu_full

    if mof in mof_isotherms:
        plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
                [(rho_hi_p-rho_lo_p)/density_units, (rho_hi_p-rho_lo_p)/density_units],
                colors.symbol(basename)+'-', label=colors.latex_mof(mof), color=colors.color(mof))

stepby = None
if 20 < ymax < 200:
    stepby = 10
elif 200 < ymax < 1000:
    stepby = 100
elif 1000 < ymax < 3000:
    stepby = 500
if stepby is not None:
    plt.yticks(list(range(0,int(ymax),stepby))
               + [round(ymax)]
               + list(range((int(ymax/stepby)+1)*stepby, int(y_max_lim), stepby)))

plt.ylim(0, 400)
plt.xlim(0, 15)
plt.yticks([])

arrowprops = dict(arrowstyle = "->")
plt.annotate('upper bound', (10, 170), xytext=(10, 250), arrowprops=arrowprops)
plt.annotate('real materials', (7, 130), xytext=(8, 80), arrowprops=arrowprops)

plt.xlabel(r'attractive energy (kJ/mol)')
plt.ylabel(r'deliverable methane')
plt.tight_layout()
plt.savefig('figs/n-vs-G-cartoon.pdf', transparent=True)
plt.savefig('figs/n-vs-G-cartoon.tiff', dpi=900)

if 'noshow' not in sys.argv:
    plt.show()
