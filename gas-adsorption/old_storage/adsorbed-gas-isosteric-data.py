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
# python3 adsorbed-gas-isosteric-data.py MOF5 absolute methane 300
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

mof = sys.argv[1]
ads_type = sys.argv[2]
basename = sys.argv[3]
my_isostere = int(sys.argv[4])

def my_interp(x0,x,y):
    return np.interp(x0, x, y)

def mu_from_n(particular_n):
    return my_interp(particular_n, n, mu)

def n_from_mu(particular_mu):
    return my_interp(particular_mu, mu, n)

def mu_from_p(particular_p):
    return my_interp(particular_p, p, mu)

data = np.loadtxt(mof + ('-%s-%s.txt' % (ads_type,basename)), skiprows=3,unpack=True,delimiter=',')

# define the temperature range for the papers data
T_label = []
Temp_str = np.loadtxt(mof + ('-%s-%s.txt' % (ads_type,basename)),dtype=str,unpack=True,delimiter=',')
for i in range(len(Temp_str)):
    if Temp_str[i][0]:
        T_label.append(int(Temp_str[i][0]))

# ---------------- Single-site Langmuir Fit (from reference) ------------------#
Lang_rho = []
Lang_p = []
Lang_T = []

Lang_T_2d, Lang_p_2d = np.meshgrid(T_label, np.arange(0, 100.1, 0.1))
Lang_rho_2d = np.zeros_like(Lang_p_2d)
Lang_rho_2d[:,:] = np.nan
axis_value = 0
# Table S3 Langmuir fitting parameters (MOF-5):
nsat1 = 30.5*22.414 # mmol/g*22.414 mL/mol --> cm3 (STP)/g
S1 = 9.2            # -R
E1 = 12.3           # kJ/mol site-specific binding energy (But this should be negative?)
R = 8.314           # J/(mol K)

for i in range(len(T_label)):
    T = T_label[i]
    # formula from paper for Single-site
    p_i = np.arange(0,100.1,0.1)
    T_i = np.full((len(p_i)),T)

    b1 = np.exp(-S1)*np.exp(E1*1000/(R*(T+273.15)))
    rho_i = nsat1*b1*p_i / (1 + b1*p_i)

    for j in range(len(p_i)):
        Lang_p.append(p_i[j])
        Lang_T.append(T_i[j])
        Lang_rho.append(rho_i[j])

    # begin interpolation of data
    Lang_p_interp = p_i
    Lang_rho_2d[:len(rho_i),axis_value] = rho_i
    axis_value += 1

    if i == 0:
        plt.plot(p_i,rho_i,'-',label=r'Langmuir', color=colors[i])
    else:
        plt.plot(p_i,rho_i,'-', color=colors[i])

plt.ylim(0,550)
plt.legend(loc = 'best')
#plt.show()
#sys.exit("Remove: Working on adding Single-Site Langmuir Model")
# ---------------- Data from the Mason and Long paper -------------------------#
#plt.figure()

isostere_T = []
isostere_p = []

all_p = []
all_T = []
all_rho = []

T_2d, p_2d = np.meshgrid(T_label, np.arange(0, 101.0, 1.0))
rho_2d = np.zeros_like(p_2d)
rho_2d[:,:] = np.nan

for i in range(0,len(data),2):
    rho_MOF = data[i+1]
    p_MOF = data[i]
    T_MOF = np.full((len(p_MOF)),T_label[i//2])
    for j in range(len(p_MOF)):
        all_p.append(p_MOF[j])
        all_rho.append(rho_MOF[j])
        all_T.append(T_MOF[j])
        #print(T_MOF[j], p_MOF[j], rho_MOF[j])

    # begin interpolation of data
    p_vals = np.arange(0,np.max(p_MOF))
    rho_interp = np.interp(p_vals,p_MOF,rho_MOF)
    rho_2d[:len(rho_interp),i//2] = rho_interp

    if np.max(rho_MOF) > my_isostere and np.min(rho_MOF) < my_isostere:
        if my_isostere < rho_MOF.max():
            # split
            where = np.argmax(rho_MOF)
            print('where is', where, rho_MOF[where], p_MOF[where])
            isostere_T.append(T_label[i//2])
            isostere_p.append(np.interp(my_isostere, rho_MOF[:where], p_MOF[:where]))
            isostere_T.append(T_label[i//2])
            isostere_p.append(np.interp(my_isostere, np.flip(rho_MOF[where:]), np.flip(p_MOF[where:])))
        else:
            isostere_T.append(T_label[i//2])
            isostere_p.append(np.interp(my_isostere, rho_MOF, p_MOF))

    plt.plot(p_MOF,rho_MOF,'.',color=colors[i//2], label = ('%s' % T_label[i//2]) + r'$^\circ C$',markersize=10)
    plt.plot(p_vals,rho_interp,linestyle='dashed',linewidth=0.5,color='black')

plt.title('%s Adsorption Data from Supplemental Info' % mof)
plt.xlabel('pressure (bar)')
plt.ylabel(r'$CH_4$ %s uptake ($\textrm{cm}^3$ (STP)/g)' % ads_type)
plt.legend(loc = 'best')
plt.savefig(mof + '-' + ads_type + '-' + basename + '-isotherm.pdf', transparent=True)

# plt.figure('rho vs T at fixed p')
# for p in np.arange(5.0, 90.0, 10):
#     rho_at_p = []
#     T_at_p = []
#     for i in range(len(T_label)):
#         #if rho_MOF:
#         isostere_p.append(np.interp(my_isostere, rho_MOF[:where], p_MOF[:where]))


# --------------- 2D Interpolation Contour Plots ------------------------------#
plt.figure()
plt.plot(isostere_T, isostere_p, 'x', label='isostere')
plt.xlabel('Temperature (K)')
plt.ylabel(r'pressure')
plt.legend(loc = 'best')

plt.figure()

def rounddown(x):
    return int(mh.floor(x / 100.0)) * 100

max_contour = rounddown(np.nanmax(rho_2d))
plt.contour(Lang_T_2d, Lang_p_2d, Lang_rho_2d, range(0,max_contour+1,50))
cbar = plt.colorbar()
cbar.set_label((r'%s uptake ($\textrm{cm}^3$ (STP)/g)' % ads_type), labelpad=10.0, y=0.6, rotation=90)
plt.title('%s isosteric Langmuir' % mof)
plt.xlabel('Temperature ($^\circ C$)')
plt.ylabel('pressure (bar)')
plt.axvline(25)
plt.savefig(mof + '-' + ads_type + '-' + basename + '-isostere.pdf', transparent=True)

plt.figure()
max_contour = rounddown(np.nanmax(rho_2d))
plt.contour(T_2d, p_2d, rho_2d, range(0,max_contour+1,50))
plt.plot(isostere_T, isostere_p, 'o', label='isostere')
cbar = plt.colorbar()
cbar.set_label((r'%s uptake ($\textrm{cm}^3$ (STP)/g)' % ads_type), labelpad=10.0, y=0.6, rotation=90)
plt.title('%s isosteric interpolation' % mof)
plt.xlabel('Temperature ($^\circ C$)')
plt.ylabel('pressure (bar)')
plt.axvline(25)
# ------------------------- Qst Calculations and Plots ------------------------#
plt.figure()

rho_range = range(50,300,5)
print(rho_range)
qst = []
Lang_qst = []
for r0 in rho_range:
    print(r0)
    contour = plt.contour(T_2d, p_2d, rho_2d, [r0])
    path = contour.collections[0].get_paths()[0]
    vert = path.vertices
    T_path = np.sort(vert[:,0])
    p_path = np.sort(vert[:,1])

    T_index = np.argwhere(T_path==25)
    delta_p = p_path[T_index+1]-p_path[T_index-1]
    delta_T = T_path[T_index+1]-T_path[T_index-1]
    dlnpdt = delta_p/delta_T/p_path[T_index]
    #dlnpdt =1/(12*4)*(np.log(p_path[T_index-2]) - 8*np.log(p_path[T_index-1]) + 8*np.log(p_path[T_index+1]) - np.log(p_path[T_index+2]))

    T2 = (25)**2
    qst_average = np.mean(T2*dlnpdt)
    print(qst_average, delta_p, delta_T, 'pressure', p_path[T_index], 'T', T_path[T_index])
    qst.append(qst_average)

    #---------------- Langmuir Part --------------------#
    Lang_contour = plt.contour(Lang_T_2d, Lang_p_2d, Lang_rho_2d, [r0])
    Lang_path = Lang_contour.collections[0].get_paths()[0]
    Lang_vert = Lang_path.vertices
    Lang_T_path = np.sort(Lang_vert[:,0])
    Lang_p_path = np.sort(Lang_vert[:,1])

    Lang_T_index = np.argwhere(Lang_T_path==25)

    Lang_delta_p = Lang_p_path[Lang_T_index+1]-Lang_p_path[Lang_T_index-1]
    Lang_delta_T = Lang_T_path[Lang_T_index+1]-Lang_T_path[Lang_T_index-1]
    Lang_dlnpdt = Lang_delta_p/Lang_delta_T/Lang_p_path[Lang_T_index]

    Lang_T2 = (25)**2
    Lang_qst_average = np.mean(Lang_T2*Lang_dlnpdt)
    Lang_qst.append(Lang_qst_average)

plt.figure()
plt.plot(rho_range,qst,'r+')
plt.plot(rho_range,Lang_qst,'o')
plt.xlabel(r'%s uptake ($\textrm{cm}^3$ (STP)/g)' % ads_type)
plt.ylabel(r'-Qst (kJ/mol)')
#plt.ylim(10,22)
plt.xlim(0,300)
plt.savefig(mof + '-' + ads_type + '-' + basename + '-Qst.pdf', transparent=True)

if 'noshow' not in sys.argv:
    plt.show()
