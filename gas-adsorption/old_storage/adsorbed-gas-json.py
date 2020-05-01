from __future__ import division, print_function
import sys, os, re, matplotlib, glob
import numpy as np
import matplotlib.pyplot as plt
import numericalunits
import json

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
mmol = numericalunits.mmol
psi = numericalunits.psi

mg = numericalunits.mg    # milligram
mL = numericalunits.mL    # milliliter
kJ = numericalunits.kJ     # kiloJoule
meter = numericalunits.m
cm = numericalunits.cm
wt = 1

mol_per_L_STP = 22.414*mol/L

# define H2 and CH4 constants
molar_mass_H2 = 2.016*gram/mol
molar_mass_CH4 = 16.04*gram/mol

rho_H2 = 8.988e-5*gram/(cm**3)
rho_CH4 = 0.000656*gram/(cm**3)

def molar_mass(json):
    masses = {
    'Hydrogen': 2.016*gram/mol,
    'H2': 2.016*gram/mol,
    'CH4': 16.04*gram/mol,
    'Methane': 16.04*gram/mol,
    }
    return masses[json['adsorbates'][0]['name']]

def crystal_density(json):
    densities = {
    'IRMOF-1': 0.61*gram/cm**3, # Furuta, Kanoya, Sakai, and Hosoe
    'CuPDC': 1.196*gram/cm**3, # Furuta, Kanoya, Sakai, and Hosoe

    #DOI:10.1021/ic201376t
    # Synonyms Basolite C300, C300, Cu-BTC, Cu3(BTC)2, HKUST-1, MOF-199
    'HKUST-1': 0.88*gram/cm**3,
    'PCN-6-Prime': 0.28*gram/cm**3, # a prime version of PCN-6?
    'PCN-6': 0.56*gram/cm**3,
    'PCN-20': 0.49*gram/cm**3,
    'MOF-399': 0.13*gram/cm**3,
    'MOF-14': 0.72*gram/cm**3,

    #DOI: 10.1126/science.1067208
    'IRMOF-1': 1.00*gram/cm**3,
    'IRMOF-3': 0.76*gram/cm**3,
    'IRMOF-14': 0.37*gram/cm**3,

    # Exceptional H2 Saturation Uptake in Microporous Metal−Organic Frameworks
    # Antek G. Wong-Foy, Adam J. Matzger*, and Omar M. Yaghi*†;
    'MOF-177': 0.477*gram/cm**3,

    'CuBTC': 0.35*gram/cm**3,


    }
    return densities[json['adsorbent']['name']]

def pressure_convert(p_unit):
    '''give each unit in units of bar'''
    if p_unit == 'bar':
        return bar
    raise ValueError

def density_units(json):
    ads_unit = json['adsorptionUnits']
    '''give each unit in units of mol/gram'''
    if ads_unit == 'cm3(STP)/g':
        return mol_per_L_STP*cm**3/gram # converts to mol/gram
    elif ads_unit == 'mmol/g':
        return 1e-3*mol/gram
    elif ads_unit == 'wt%':
        return (1/molar_mass(json))
    elif ads_unit == 'mg/g':
        return (1e-3/molar_mass(json))
    else:
        raise ValueError

    # 'cm^3/g': cm**3/gram, need to know the temperature
    # 'cm^3/cm^3': cm**3/cm**3, need to know the density of the material

def extract_units(json):
    T_units = Kelvin
    P_units = json['pressureUnits']

    # return the Temperature, Pressure, and Total Adsorption units
    P_units = pressure_convert(P_units)
    ads_units = density_units(json)


    return T_units, P_units, ads_units

def stored_data(json, adsorbate):
    P = []
    ads = []
    T = json['temperature']
    for pressure in json['isotherm_data']:
        P.append(pressure['pressure'])

    for ads_total in json['isotherm_data']:
        ads_add = ads_total['species_data']
        ads_add = list(map(lambda x: x["adsorption"], ads_add))
        ads.append(ads_add[0])

    # return the Temperature, Pressure, and Total Adsorption
    T_units, P_units, ads_units = extract_units(json)
    print('T_units', T_units/Kelvin, 'K')
    print('P_units', P_units/bar, 'bar')
    print('ads_units', ads_units, 'mol/gram')
    return T*T_units, np.array(P)*P_units, np.array(ads)*ads_units

def read_experimental_data(adsorbate): # the adsorbate is neccessary in order to do unit conversion
    ads_dict = {}
    gas = adsorbate
    if adsorbate == 'Hydrogen':
        adsorbate = 'H2'
    if adsorbate == 'Methane':
        adsorbate = 'CH4'

    for json_name in glob.glob('*.json'):
        with open('%s' % json_name) as json_file:
            #try:
                json_data = json.load(json_file)
                if json_data['adsorbates'][0]['name'] == gas:
                    print('The filename is', json_name)
                    try:
                        crysden = crystal_density(json_data)
                    except:
                        print('unknown density for', json_data['adsorbent']['name'])
                        continue
                    T, p, ads = stored_data(json_data, adsorbate)
                    ads_dict['%s' % json_name.split('.')[0]] = {
                    'temperature': T,
                    'pressure': p,
                    'ads': ads,
                    'rho': ads*crysden, # in mol/volume
                    'crystal_density': crysden, # in mass/volume
                    }
            #except:
            #    print('unable to handle', json_name)
    return ads_dict
