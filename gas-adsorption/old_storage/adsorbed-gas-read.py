from __future__ import division, print_function
import sys, os, re, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import numericalunits

matplotlib.rcParams['text.usetex'] = True

import colors

from matplotlib.font_manager import FontProperties

small_font = FontProperties()
small_font.set_size('small')

# import python json data extract package
ads_json = __import__('adsorbed-gas-json')

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

basename = sys.argv[1]
basename_at_dif_T = sys.argv[2]
p_empty = float(sys.argv[3])*bar
p_full = float(sys.argv[4])*bar


molecular_weight = 2.016*gram/mol
if 'methane' in basename:
    molecular_weight = 16.04*gram/mol

kB = 8.314*Joule/mol/Kelvin
data = [np.loadtxt('%s.csv' % basename, skiprows=1),
             np.loadtxt('%s.csv' % basename_at_dif_T, skiprows=1)]

T = data[0][:,0]*Kelvin            # Temperature
#print(np.asarray(T).shape)
p = data[0][:,1]*bar               # Pressure
rho = data[0][:,2]*mol/L           # Density
V = data[0][:,3]*L/mol             # Volume
U = data[0][:,4]*kJ/mol            # Internal Energy
H = data[0][:,5]*kJ/mol            # Enthalpy
S = data[0][:,6]*Joule/mol/Kelvin  # Entropy

G = H - T*S
F = U - T*S

mu = G # per mol
n = 1/V
kT = T[0]*kB

T2 = data[1][:,0]*Kelvin            # Temperature
p2 = data[1][:,1]*bar               # Pressure
rho2 = data[1][:,2]*mol/L           # Density
V2 = data[1][:,3]*L/mol             # Volume
U2 = data[1][:,4]*kJ/mol            # Internal Energy
H2 = data[1][:,5]*kJ/mol            # Enthalpy
S2 = data[1][:,6]*Joule/mol/Kelvin  # Entropy

G2 = H2 - T2*S2
F2 = U2 - T2*S2

mu2 = G2 # per mol
n2 = 1/V2
kT2 = T2[0]*kB

max_Gads = 60 # Just pick a Big number

def my_interp(x0,x,y):
    return np.interp(x0, x, y)

def mu_from_n(particular_n):
    return my_interp(particular_n, n, mu)

def n_from_mu(particular_mu):
    return my_interp(particular_mu, mu, n)

def mu_from_p(particular_p):
    return my_interp(particular_p, p, mu)

mu_empty = mu_from_p(p_empty)
mu_full = mu_from_p(p_full)

def n_from_mu2(particular_mu):
    return my_interp(particular_mu, mu2, n2)

def mu_from_p2(particular_p):
    return my_interp(particular_p, p2, mu2)

mu_empty2 = mu_from_p2(p_empty)
mu_full2 = mu_from_p2(p_full)

if 'methane' in basename:
    density_units = 0.044135*mol/L # density at STP in mol/L
    density_unit_name = 'cm$^3$ CH$_4$ STP/cm$^3$'
    gas = 'Methane'
else:
    #density_units = kg/L/molecular_weight
    #density_unit_name = 'kg H$_2$/L'
    density_units = gram/L/molecular_weight
    density_unit_name = 'g H$_2$/L'
    gas = 'Hydrogen'

Gads = np.linspace(0,max_Gads,1000)*kJ/mol

experimental_data = {}
# DATA is from TABLE 3 DOI: 10.1021/acsenergylett.8b00154
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
    'COF102': 0.41*gram/cm**3, # DOI: 10.1039/B805473H
}

H2_298_298_Data = {
    'HKUST1': {
        'rho_c': 6.08*mmol/gram,
        'rho_h': 0.36*mmol/gram,
    },
    'NOTT112': {
        'rho_c': 8.83*mmol/gram,
        'rho_h': 0.47*mmol/gram,
    },
    'NU125': {
        'rho_c': 8.16*mmol/gram,
        'rho_h': 0.51*mmol/gram,
    },
    'rhtMOF7': {
        'rho_c': 5.72*mmol/gram,
        'rho_h': 0.34*mmol/gram,
    },
    'CuMOF74': {
        'rho_c': 3.78*mmol/gram,
        'rho_h': 0.22*mmol/gram,
    },
    'PCN250': {
        'rho_c': 5.75*mmol/gram,
        'rho_h': 0.34*mmol/gram,
    },
    'NU1000': {
        'rho_c': 8.67*mmol/gram,
        'rho_h': 0.51*mmol/gram,
    },
    'UiO67': {
        'rho_c': 5.74*mmol/gram,
        'rho_h': 0.33*mmol/gram,
    },
    'UiO68Ant': {
        'rho_c': 7.31*mmol/gram,
        'rho_h': 0.38*mmol/gram,
    },
    'CYCU3Al': {
        'rho_c': 8.91*mmol/gram,
        'rho_h': 0.49*mmol/gram,
    },
    'Zn2bdc2dabco': {
        'rho_c': 5.29*mmol/gram,
        'rho_h': 0.31*mmol/gram,
    },
}

H2_77_160_Data = {
    'HKUST1': {
        'qh': 6.3*kJ/mol,
        'qc': 4*kJ/mol,
        'rho_c': 28.04*mmol/gram,
        'rho_h': 2.23*mmol/gram,
    },
    'NOTT112': {
        'qh': 4.8*kJ/mol,
        'qc': 47.86*kJ/mol,
        'rho_c': 26*mmol/gram,
        'rho_h': 2.54*mmol/gram,
    },
    'NU125': {
        'qh': 4.8*kJ/mol,
        'qc': 3.5*kJ/mol,
        'rho_c': 44.67*mmol/gram,
        'rho_h': 2.62*mmol/gram,
    },
    'rhtMOF7': {
        'qh': 5.9*kJ/mol,
        'qc': 3.2*kJ/mol,
        'rho_c': 25.68*mmol/gram,
        'rho_h': 2.13*mmol/gram,
    },
    'CuMOF74': {
        'qh': 5.8*kJ/mol,
        'qc': 4*kJ/mol,
        'rho_c': 16.26*mmol/gram,
        'rho_h': 1.51*mmol/gram,
    },
    'PCN250': {
        'qh': 6.2*kJ/mol,
        'qc': 4.4*kJ/mol,
        'rho_c': 28.35*mmol/gram,
        'rho_h': 2.43*mmol/gram,
    },
    'NU1000': {
        'qh': 4.7*kJ/mol,
        'qc': 1.6*kJ/mol,
        'rho_c': 43.38*mmol/gram,
        'rho_h': 2.115*mmol/gram,
    },
    'UiO67': {
        'qh': 5.5*kJ/mol,
        'qc': 3*kJ/mol,
        'rho_c': 31.14*mmol/gram,
        'rho_h': 1.85*mmol/gram,
    },
    'UiO68Ant': {
        'qh': 5.5*kJ/mol,
        'qc': 3*kJ/mol,
        'rho_c': 41*mmol/gram,
        'rho_h': 2.34*mmol/gram,
    },
    'CYCU3Al': {
        'qh': 4.4*kJ/mol,
        'qc': 2*kJ/mol,
        'rho_c': 44.92*mmol/gram,
        'rho_h': 1.90*mmol/gram,
    },
    'Zn2bdc2dabco': {
        'qh': 4.9*kJ/mol,
        'qc': 4.5*kJ/mol,
        'rho_c': 25.72*mmol/gram,
        'rho_h': 1.65*mmol/gram,
    },
}

H2_77_77_Data = { # from Table S4. In the supplemental.
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
    'COF102': {
        'qh': 3.5*kJ/mol,  #-> 4.7*kJ/mol,
        'qc': 3.5*kJ/mol,  # only have low pressure so don't know?
        'rho_c': 62*mg/gram/molecular_weight, #-> 25.72*mmol/gram,
        'rho_h': 43*mg/gram/molecular_weight, #-> 17.93*mmol/gram,
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
    # include the COF-102 data from 10.1021/ja9015765
    'COF102': 12*kJ/mol,
}


# Note: I took cm^3/g and converted to L/mol using molecular weight
# DOI:10.1039/C3SC52633J
COF102_Crys_Dens = 0.17*gram/cm**3
COF_density_units = mg/gram * COF102_Crys_Dens / molecular_weight
print('The COF DENSITY IS',COF_density_units)
def linear_interpolate(x, x1, y1, x2, y2):
    return ((x2-x)*y1 + (x-x1)*y2)/(x2-x1)
if basename == 'methane-298' and basename_at_dif_T == 'methane-298':
    experimental_data = {
        'HKUST1':
        (  linear_interpolate(65.5, 65.05, 263.01, 70.01, 266.18)*density_units
         - linear_interpolate(5.8 , 5.57 , 79.60 ,  7.71, 101.94)*density_units,
           17*kJ/mol),
        'PCN14':
        (  linear_interpolate(65.5, 64.89, 239.42, 69.91, 243.67)*density_units
         - linear_interpolate(5.8 , 5.58, 81.23, 7.56, 100.16)*density_units,
           17.6*kJ/mol),
        'Ni2dobdc':
        (  linear_interpolate(65.5, 65.15, 259.64, 70.06, 260.58)*density_units
         - linear_interpolate(5.8 , 5.79, 126.15, 7.73 , 144.71)*density_units,
           20.7*kJ/mol),
        'Co2dobdc':
        (  linear_interpolate(65.5, 64.89, 248.64, 69.63, 252.11)*density_units
         - linear_interpolate(5.8 , 5.62, 117.78, 7.64, 137.04)*density_units,
           19.5*kJ/mol),
        'Mg2dobdc':
        (  linear_interpolate(65.5, 64.90, 230.28, 69.93, 233.71)*density_units
         - linear_interpolate(5.8 , 5.64, 93.50, 7.60, 111.58)*density_units,
           18.5*kJ/mol),
        'MOF5':
        (  linear_interpolate(65.5, 64.91, 215.14, 69.90, 222.03)*density_units
         - linear_interpolate(5.8 , 5.58, 28.97, 7.57, 39.00)*density_units,
           12.2*kJ/mol),
        'AX21':
        (  linear_interpolate(65.5, 64.87, 202.80, 69.65, 209.63)*density_units
         - linear_interpolate(5.8 , 5.55, 52.19, 7.53, 64.41)*density_units,
           14*kJ/mol),
        # include the COF-102 data from 10.1021/ja9015765
        # I used WebPlotDigitizer and extracted the data now interpolate!
        'COF102':
        (  linear_interpolate(65.5, 59.49, 231.87, 66.03, 237.75)*COF_density_units
         - linear_interpolate(5.8 , 4.07, 27.65, 6.73, 43.52)*COF_density_units,
           12*kJ/mol),
    }

plt.figure('deliverable capacity', figsize=(5,3.75))

D = n_from_mu(mu_full+Gads)-n_from_mu2(mu_empty2+Gads)

y = D/density_units
ymax = max(y)
xpos = np.where(y == ymax)[0]
x = Gads/(kJ/mol)
xmax = x[xpos][0]

XX = Gads/(kJ/mol)
YY = (n_from_mu(mu_full+Gads)-n_from_mu(mu_empty+Gads))/density_units

max_index = np.argwhere(YY == np.amax(YY))[0][0]
derivative = np.diff(YY)[max_index:]/np.diff(XX)[max_index:]
minval = np.nanmin(derivative[np.isfinite(derivative)])
new_index = np.argwhere(derivative != 0)[-10][0]

biggest_x_reasonable = Gads[0:new_index].max()/(kJ/mol)
plt.xlim(0,biggest_x_reasonable)
plt.ylim(0,np.max(n_from_mu(mu_full+Gads))/density_units)

x_max_lim = min(3*xmax, biggest_x_reasonable)
plt.xlim(0,x_max_lim)
if basename == 'methane-298':
    plt.ylim(0,1.5*ymax)
else:
    plt.ylim(0,1.5*ymax)

to_be_legended = []
legend_labels = []
# plt.title(r'%s $p_{\mathrm{full}} = $ %g, $p_{\mathrm{empty}} = $ %g' % (gas,p_full/bar,p_empty/bar))
to_be_legended.append(
    plt.plot(Gads/(kJ/mol),
             n_from_mu(mu_full+Gads)/density_units, ':')[0])
legend_labels.append(r'$\rho(\mu_{\mathrm{full}} + \Delta g_{st}$)')
# plt.plot(Gads/(kJ/mol),
#          n_from_mu(mu_empty+Gads)/density_units, ':', label = r'$\rho_{L}$')
# D = n_from_mu(mu_full+Gads)-n_from_mu(mu_empty+Gads)
# plt.plot(Gads/(kJ/mol), D/density_units, label = r'$D$')

to_be_legended.append(
    plt.plot(Gads/(kJ/mol),
             n_from_mu2(mu_empty2+Gads)/density_units, ':')[0])
legend_labels.append(r'$\rho(\mu_{\mathrm{empty}} + \Delta g_{st}$)')

to_be_legended.append(plt.plot(Gads/(kJ/mol), D/density_units)[0])
legend_labels.append(r'$\rho_D(\Delta g_{st})$')

n_goals = []
n_arpa_methane_goal = 315*0.044135*mol/L # 315* density at STP
if 'methane' in basename:
    print(n_arpa_methane_goal/(mol/L))
    n_goals.append((n_arpa_methane_goal, 'ARPA-E', '#999999', '-'))
else:
    n_goals.append((0.05*kg/L/molecular_weight, 'DOE ULTIMATE','#999999', '-'))
    n_goals.append((0.04*kg/L/molecular_weight, 'DOE 2025',    '#888888', '-'))
    n_goals.append((0.03*kg/L/molecular_weight, 'DOE 2020',    '#777777', '-'))

bbox = dict(boxstyle="round", fc="0.8")
arrowprops = dict(arrowstyle = "-")
# plt.plot(xmax,ymax,'.',mfc='0.8',color = 'k',markersize=4)
# plt.annotate('(%.3f, %.3f)'%(xmax, ymax),
#             (xmax, ymax), xytext=(10,10),
#             xycoords='data',
#             textcoords='offset points',
#             bbox=bbox, arrowprops=arrowprops)
# max_label = '%.2g-' % ymax
# if ymax >= 10:
#     max_label = r'%.0f -' % round(ymax)
# plt.text(0, ymax, max_label,
#          horizontalalignment='right',
#          va='center',
#          # bbox=dict(facecolor='white', edgecolor='white', alpha=1, pad=0)
# )
plt.axhline(ymax, color='g', linestyle='--', linewidth=0.5)


x_goal_label = 0.03*x_max_lim
for n_goal, label, color, style in n_goals:
    plt.text(x_goal_label, n_goal/density_units+0.02*ymax, label, color=color,
             bbox=dict(facecolor='white', edgecolor='white', alpha=1, pad=0))
    line = plt.axhline(n_goal/density_units, color=color, linestyle=style,
                       linewidth=0.5)
    # to_be_legended.append(line)
    # legend_labels.append(label)

plt.text(x_goal_label, 1.2*ymax,
         '%s\n$p_{\mathrm{full}} = $ %g bar\n$p_{\mathrm{empty}} = $ %g bar\n$T=%g$ K'
         % (gas,
            p_full/bar,
            p_empty/bar,
            T[0]/Kelvin,
         ),
         color='k')

first_legend = plt.legend(to_be_legended, legend_labels, loc='upper right')
plt.gca().add_artist(first_legend)

for mof in colors.order(experimental_data):
    if basename == 'methane-298':
        try:
            data = np.loadtxt(mof +'-298.txt', skiprows=1,unpack=True)
            T_MOF = data[0][0]*Kelvin
            p_MOF = data[1]*bar
            print(mof)
            rho_MOF = data[2]*density_units

            rho_lo_p = my_interp(p_empty, p_MOF, rho_MOF)
            delta_G_lo_p = mu_from_n(rho_lo_p) - mu_from_p(p_empty)

            rho_hi_p = my_interp(p_full, p_MOF, rho_MOF)
            delta_G_hi_p = mu_from_n(rho_hi_p) - mu_from_p(p_full)
            plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
                    2*[experimental_data[mof][0]/density_units],
                    '.-', label=colors.latex_mof(mof), color=colors.color(mof))
            # plt.plot(experimental_data[mof][1]/(kJ/mol),
            #          experimental_data[mof][0]/density_units, '.')
            rho_MOF_empty = my_interp(p_empty, p_MOF, rho_MOF)
            print('    rho empty', mof, rho_MOF_empty/density_units)
            print('delta H', mof,
                  (my_interp(rho_MOF_empty, n, H) - my_interp(p_empty, p, H))/(kJ/mol))
        except:
            pass
    # else:
    #     plt.plot(experimental_data[mof][1]/(kJ/mol),
    #              experimental_data[mof][0]/density_units,
    #              'x', label=mof, color=colors.color(mof))

for mof in colors.order(H2_77_160_Data): # For each MOF
    if basename == 'hydrogen-77' and basename_at_dif_T == 'hydrogen-160':
        rho_lo_p = H2_77_160_Data[mof]['rho_h']*H2_crystal_density[mof]
        delta_G_lo_p = my_interp(rho_lo_p, n2, mu2) - mu_empty2
        # convert units: mmol/gram*molecular_weight/crystal_density,
        rho_hi_p = H2_77_160_Data[mof]['rho_c']*H2_crystal_density[mof]
        delta_G_hi_p = my_interp(rho_hi_p, n, mu) - mu_full

        plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
                [(rho_hi_p-rho_lo_p)/density_units, (rho_hi_p-rho_lo_p)/density_units],
                'x-', label=colors.latex_mof(mof), color=colors.color(mof))

for mof in colors.order(H2_77_77_Data): # For each MOF
    if basename == 'hydrogen-77' and basename_at_dif_T == 'hydrogen-77':
        rho_lo_p = H2_77_77_Data[mof]['rho_h']*H2_crystal_density[mof]
        delta_G_lo_p = my_interp(rho_lo_p, n2, mu2) - mu_empty2
        # convert units: mmol/gram*molecular_weight/crystal_density,
        rho_hi_p = H2_77_77_Data[mof]['rho_c']*H2_crystal_density[mof]
        delta_G_hi_p = my_interp(rho_hi_p, n, mu) - mu_full

        plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
                [(rho_hi_p-rho_lo_p)/density_units, (rho_hi_p-rho_lo_p)/density_units],
                'x-', label=colors.latex_mof(mof), color=colors.color(mof))

for mof in colors.order(H2_298_298_Data): # For each MOF
    if basename == 'hydrogen-298' and basename_at_dif_T == 'hydrogen-298':
        rho_lo_p = H2_298_298_Data[mof]['rho_h']*H2_crystal_density[mof]
        delta_G_lo_p = my_interp(rho_lo_p, n2, mu2) - mu_empty2
        # convert units: mmol/gram*molecular_weight/crystal_density,
        rho_hi_p = H2_298_298_Data[mof]['rho_c']*H2_crystal_density[mof]
        delta_G_hi_p = my_interp(rho_hi_p, n, mu) - mu_full

        plt.plot([delta_G_lo_p/(kJ/mol), delta_G_hi_p/(kJ/mol)],
                [(rho_hi_p-rho_lo_p)/density_units, (rho_hi_p-rho_lo_p)/density_units],
                'x-', label=colors.latex_mof(mof), color=colors.color(mof))

print('delta H', (my_interp(mu_empty + xmax, mu, H) - my_interp(p_empty, p, H))/(kJ/mol))
print('delta H', (my_interp(mu_full + xmax, mu, H) - my_interp(p_full, p, H))/(kJ/mol))

plt.xlabel(r'$\left|\Delta g_{st}\right|$ (kJ/mol)')
plt.ylabel(r'$\rho_D$ (%s)' % density_unit_name)
if ymax < 300:
    plt.yticks([0,10,20,30,40,50,round(ymax),60,70,80])
else:
    plt.yticks([0,100,200,300,round(ymax),400,500,600],
               labels=[0,100,200,300,'%.0f' % round(ymax),400,500,600])
plt.legend(loc = 'lower right', framealpha=1.0, prop=small_font)
plt.tight_layout()
plt.savefig(basename + '-' + basename_at_dif_T + '-n-vs-G.pdf', transparent=True)


plt.figure()
for n_goal, label, color, style in n_goals:
    plt.plot(Gads/(kJ/mol), n_goal/D, label = label, color=color, linestyle=style)
plt.xlabel(r'$\Phi$ (kJ/mol)')
plt.ylabel(r'$\epsilon_{\min}$')
plt.xlim(0,Gads[0:new_index].max()/(kJ/mol))
plt.ylim(0,1)
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig(basename + '-' + basename_at_dif_T + '-void.pdf', transparent=True)

#plt.figure()

#def best_n_deliv(p_empty, p_full):
    #mu_empty = mu_from_p(p_empty)
    #mu_full = mu_from_p(p_full)
    #Gads = np.linspace(0*kJ/mol,50*kJ/mol,10000)
    #n_deliv = n_from_mu(mu_full+Gads)-n_from_mu(mu_empty+Gads)
    #return np.amax(n_deliv)

#FULL, EMPTY = np.meshgrid(np.linspace(1, 100, 100)*bar,
                          #np.linspace(1, 100, 100)*bar)
#best = np.zeros_like(EMPTY)
#for i in range(len(EMPTY)):
    #for j in range(len(EMPTY)):
        #best[i,j] = best_n_deliv(EMPTY[i,j], FULL[i,j])
#plt.contour(EMPTY, FULL, best, 100, linewidths=1)
#plt.contourf(EMPTY, FULL, best, 100)

#cbar = plt.colorbar()
#cbar.set_label('n_deliv', labelpad=-40, y=1.05, rotation=0)

#plt.xlabel(r'$p_{empty}$ (bar)')
#plt.ylabel(r'$p_{full}$ (bar)')

#plt.savefig(basename + '-contour-plot.pdf', transparent=True)
# plt.figure()
# plt.title('Ratio of Heat of Adsorptions')
# plt.plot(p/(bar),p/kT/density_units * (rho - rho[0]*np.exp((mu-mu[0])/kT))/density_units)
# plt.xlim(0, 100)
# plt.ylim(0, 3000)
# plt.xlabel(r'pressure (bar)')
# plt.ylabel(r'$q_{st}$' + '/' r'$q_{st}^{expt}$')
# plt.legend(loc='best')

fig_rho_mu = plt.figure('rho vs mu', figsize=(5,3.75))
if 'methane' in basename:
    expdata = ads_json.read_experimental_data('Methane')
else:
    expdata = ads_json.read_experimental_data('Hydrogen')
mu_ax = fig_rho_mu.add_subplot(111)
p_ax = mu_ax.twiny()

full_color = '#999999' # '#333333'
empty_color = '#999999'
if kT == kT2:
    #plt.figure('rho vs mu')
    rho_max = np.max(n_from_mu(mu_full+Gads))
    # plt.title(r'%s density vs mu' % (basename))
    mu_ax.plot(mu/(kJ/mol),rho/density_units, label=r'bulk methane')
    mu_ax.plot(mu/(kJ/mol),rho[0]*np.exp((mu-mu[0])/kT)/density_units, '-',
             label='ideal gas')

    if basename == 'methane-298':
        for mof in experimental_data: # For each MOF
            try:
                data = np.loadtxt(mof +'-298.txt', skiprows=1,unpack=True)
                T_MOF = data[0][0]*Kelvin
                p_MOF = data[1]*bar
                rho_MOF = data[2]*density_units
                PHI = 0 # np.amax((mu_from_n(rho_MOF) - mu_from_p(p_MOF))/(kJ/mol))
                mu_ax.plot(mu_from_p(p_MOF)/(kJ/mol) + PHI,
                         rho_MOF/density_units, '--', label=colors.latex_mof(mof))
            except OSError: # Raise an OSError if file does not exist?
                print('OSError: The file does not exist continuing...')
                continue

    for mof in sorted(expdata.keys()):
        print(mof)
        T_MOF = expdata[mof]['temperature']
        if abs(T_MOF - kT/kB)/Kelvin > 0.5:
            print('Wrong temperature:', T_MOF/Kelvin,'vs', kT/Kelvin/kB)
            continue
        p_MOF = expdata[mof]['pressure']
        rho_MOF = expdata[mof]['rho']

        print('   T', T_MOF/Kelvin)
        print('   p', p_MOF/bar)
        print('   n', rho_MOF)
        PHI = 0 # np.amax((mu_from_n(rho_MOF) - mu_from_p(p_MOF))/(kJ/mol))
        # plt.plot(mu_from_p(p_MOF)/(kJ/mol) + PHI,
        #         rho_MOF/density_units, '-', label=mof)
        print('  mu_from_p',mu_from_p(p_MOF)/(kJ/mol) + PHI)
        print(' rho_Mof', rho_MOF)

    mu_ax.axvline(mu_from_p(p_full)/(kJ/mol), linestyle=':', color=full_color)
    mu_ax.axvline(mu_from_p(p_empty)/(kJ/mol), linestyle=':', color=empty_color)
    mu_ax.set_xlabel(r'$\mu$ (kJ/mol)')
    mu_ax.set_ylabel(r'$\rho$ (%s)' % density_unit_name)
    mu_ax.set_ylim(0, 600)
    mu_ax.set_xlim(-17, 6)
    p_ax.set_xlim(-17, 6)
    pressure_ticks = np.array([2,
                               5.8,
                               20,
                               65,
                               200,
                               1000,
                               2000,
    ])
    p_ax.set_xlabel(r'$p$ (bar)')
    p_ax.set_xticks(mu_from_p(pressure_ticks*bar)/(kJ/mol))
    p_ax.set_xticklabels(['%g' % x for x in pressure_ticks])
    mu_ax.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(basename + '-rho-mu.pdf', transparent=True)

    plt.figure()
    # plt.title(r'%s density vs pressure' % (basename))
    plt.plot(p/(bar),rho/density_units)
    plt.plot(p/bar,p/kT/density_units,':')
    plt.ylim(0, rho_max/density_units)
    plt.xlim(0, 1000)
    plt.axvline(p_full/bar, linestyle=':', color=full_color)
    plt.axvline(p_empty/bar, linestyle=':', color=empty_color)
    plt.xlabel(r'pressure (bar)')
    plt.ylabel(r'$\rho$ (%s)' % density_unit_name)
    plt.tight_layout()
    plt.savefig(basename + '-rho-pressure.pdf', transparent=True)

    if basename == 'methane-298':
        plt.figure('Gibbs vs rho')
        for mof in experimental_data: # For each MOF
            try:
                data = np.loadtxt(mof +'-298.txt', skiprows=1,unpack=True)
                T_MOF = data[0][0]*Kelvin
                p_MOF = data[1]*bar
                rho_MOF = data[2]*density_units
                plt.plot(rho_MOF/density_units,
                         (mu_from_n(rho_MOF) - mu_from_p(p_MOF))/(kJ/mol),
                         '-', label=mof)
            except OSError: # Raise an OSError if file does not exist?
                print('OSError: The file does not exist continuing...')
                continue
        plt.xlabel(r'$\rho$ (%s)' % density_unit_name)
        plt.ylabel(r'$\Delta G/N$ (kJ/mol)')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(basename + '-experimental-gibbs.pdf', transparent=True)

if 'noshow' not in sys.argv:
    plt.show()
