from __future__ import division, print_function
import sys, os, re, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import numericalunits
import math as mh
import collections
# Data from the work:  Benchmark Study of Hydrogen Storage in
# Metal−Organic Frameworks under Temperature and Pressure Swing Conditions

# HOW TO RUN:
# python3 new-adsorbed-gas-isosteric.py MOF5 absolute methane
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')

# TO DO:

# 1. Add implementation that does not assume an ideal gas for methane.
#    (Requires NIST data.)

# 3. Attempt a Gibbs free energy analysis.

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

mof = sys.argv[1]
ads_type = sys.argv[2]
basename = sys.argv[3]

def my_interp(x0,x,y):
    return np.interp(x0, x, y)

data = np.loadtxt(mof + ('-%s-%s.txt' % (ads_type,basename)), skiprows=2,unpack=True,delimiter=',')

# define the temperature range for the papers data
T_label = []
Temp_str = np.loadtxt(mof + ('-%s-%s.txt' % (ads_type,basename)),dtype=str,unpack=True,delimiter=',')
print(mof + ('-%s-%s.txt' % (ads_type,basename)))
for i in range(len(Temp_str)):
    if Temp_str[i][0]:
        T_label.append(int(Temp_str[i][0]))
T_raw = np.array(T_label)

#print(data)

pressure = data[::2,:]
rho = data[1::2,:]
T = np.zeros_like(pressure)
for i in range(len(T_raw)):
    T[i,:] = T_raw[i]

for i in range(T.shape[0]):
    plt.plot(pressure[i,:], rho[i,:], '.-', label='{}$^\circ C$'.format(T[i,0]))
plt.legend(loc='best')
plt.xlabel('$p$')
plt.ylabel(r'$\rho$')

plt.figure('rho vs T at fixed p')

for p in np.arange(5.0, 90.0, 20):
    rho_at_p = []
    T_at_p = []
    for i in range(T.shape[0]):
        if pressure[i,:].max() > p:
            rho_at_p.append(np.interp(p, pressure[i,:], rho[i,:]))
            T_at_p.append(T[i,0])
    plt.plot(T_at_p, rho_at_p, 'o-', label='{} bar'.format(p))
    print('pressure', p)
    print('    T', T_at_p)
    print('  rho', rho_at_p)
plt.legend(loc='best')
plt.xlabel('$T$')
plt.ylabel(r'$\rho$')

my_pressures = np.arange(2., 95, 0.1)
my_drho_dT = np.zeros_like(my_pressures)
my_drho_dp = np.zeros_like(my_pressures)
my_rho = np.zeros_like(my_pressures)
for j in range(len(my_pressures)):
    p = my_pressures[j]
    dp = min(p/2, 2.5)

    rho_at_p = []
    T_at_p = []
    for i in range(3): # assume 25C is the second temperature
        if pressure[i,:].max() > p:
            rho_at_p.append(np.interp(p, pressure[i,:], rho[i,:]))
            T_at_p.append(T[i,0])

    # The following takes the derivative drho_dp of a quadratic
    # interpolation through the three nearest pressure values at this
    # temperature.
    idx = (np.abs(pressure[i,:] - p)).argmin()
    if idx == 0:
        idx = 1
    if idx == pressure.shape[1]-1:
        idx = pressure.shape[1]-2
    print('idx is',idx)
    print('p', p)
    print('value is', rho[i,idx-1]*(2*p - pressure[i,idx] - pressure[i,idx+1])/(pressure[i,idx-1]-pressure[i,idx])/(pressure[i,idx-1]-pressure[i,idx+1]))
    my_drho_dp[j] = rho[i,idx-1]*(2*p - pressure[i,idx] - pressure[i,idx+1])/(pressure[i,idx-1]-pressure[i,idx])/(pressure[i,idx-1]-pressure[i,idx+1])
    my_drho_dp[j] += rho[i,idx]*(2*p - pressure[i,idx-1] - pressure[i,idx+1])/(pressure[i,idx]-pressure[i,idx-1])/(pressure[i,idx]-pressure[i,idx+1])
    my_drho_dp[j] += rho[i,idx+1]*(2*p - pressure[i,idx] - pressure[i,idx-1])/(pressure[i,idx+1]-pressure[i,idx])/(pressure[i,idx+1]-pressure[i,idx-1])

    my_rho[j] = rho_at_p[1]
    dT_plus = T_at_p[2] - T_at_p[1]
    dT_minus = T_at_p[1] - T_at_p[0]
    shifted_rho_plus = rho_at_p[2] - rho_at_p[1]
    shifted_rho_minus = rho_at_p[0] - rho_at_p[1]
    my_drho_dT[j] = ((shifted_rho_plus*dT_minus/dT_plus - shifted_rho_minus*dT_plus/dT_minus)
                     /(dT_plus + dT_minus))

plt.figure('drho/dT at 25C as a function of rho')
plt.plot(my_rho, my_drho_dT, '-')
plt.legend(loc='best')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\left(\frac{\partial\rho}{\partial T}\right)_p$')

plt.figure('drho/dp at 25C as a function of rho')
plt.plot(my_rho, my_drho_dp, '-')
plt.legend(loc='best')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\left(\frac{\partial\rho}{\partial p}\right)_T$')

plt.figure('dp/dT at 25C as a function of rho')
plt.plot(my_rho, -my_drho_dT/my_drho_dp, '-')
plt.legend(loc='best')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\left(\frac{\partial p}{\partial T}\right)_\rho$')

plt.figure('Qst at 25C as a function of rho')
#Langmuir Fit
Langmuir = np.full_like(my_rho,12.3)

plt.plot(my_rho, -25**2*my_drho_dT/my_drho_dp/my_pressures, '-')
plt.plot(my_rho, Langmuir, '--',label=r'Langmuir')
plt.legend(loc='best')
plt.xlabel(r'$\rho$')
plt.ylim(0,22)
plt.ylabel(r'$-\left(\frac{\partial p}{\partial T}\right)_\rho$')
plt.savefig(mof + '-' + ads_type + '-' + basename + '-Qst.pdf', transparent=True)


if 'noshow' not in sys.argv:
    plt.show()
