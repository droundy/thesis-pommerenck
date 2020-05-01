from __future__ import division, print_function
import sys, os, re, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import numericalunits

import colors

matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')

# All data comes from DOI: 10.1021/acsenergylett.8b00154
v_coefs = {
    'HKUST1': {
        'a': [-784.62, 11.334, 0.64936, -0.066476, 0.003472, -0.00014164, 3.71e-6],
        'b': [5.6102, -0.039763, 0.0071242, -0.00095632, 6.27e-5, -1.06e-6, -1.88e-8]
    },
    'NOTT112': {
        'a': [-613.79, 13.013, -0.21045, -0.0071214, 0.00048134, -1.0052e-5, 2.4781e-7, -5.3734e-9],
        'b': [4.431, -0.049712, 0.0040162, -0.00019197, 5.7258e-6, -6.8488e-8, -2.1531e-9, 7.1793e-11]
    },
    'NU125': {
        'a': [-593.6, -0.029702, 0.50094, -0.0088905, -0.00056209, 1.9733e-5, -1.1469e-7, 1.5171e-9],
        'b': [4.3146, 0.0030633, 0.0036414, -0.00044774, 2.0118e-5, -2.5e-7, -2.7021e-9, 2.3075e-11]
    },
    'rhtMOF7': {
        'a': [-712.95, -1.2764, 1.8356, -0.054302, -0.0026518, 0.00013347],
        'b': [5.1841, 0.055013, -0.00070978, -0.00092537, 8.6874e-5, -2.2699e-6]
    },
    'CuMOF74': {
        'a': [-703.34, -4.9553, 2.4354, -0.01936, -0.0034717, -8.9968e-5],
        'b': [5.513, 0.05651, 0.015549, -0.0036854, 0.00022566, -1.1167e-6]
    },
    'PCN250': {
        'a': [-792.52, 27.244, 0.092688, -0.49499, 0.048062, -0.001344],
        'b': [5.4822, -0.056556, -0.021315, 0.0074926, -0.00064637, 1.7667e-5]
    },
}




IsothermT = [77,298]
x = np.arange(0.001,50,0.001)

for mof in sorted(v_coefs):
    a_coef = 0
    b_coef = 0
    l_max = 0
    print(mof)
    for i in range(0,len(v_coefs[mof]['a'])):
        a_coef += v_coefs[mof]['a'][i]*x**i
    for i in range(0,len(v_coefs[mof]['b'])):
        b_coef += v_coefs[mof]['b'][i]*x**i
    for T in IsothermT:
        p = np.exp(np.log(x) + 1/T * a_coef + b_coef)

        index = np.where(p < 100)[0]
        #print(index)
        max_index = np.max(np.where(p < 100))
        for j in range(0,max_index-1):
            if (index[j+1] - index[j] > 1):
                l_max = j
                break
            else:
                l_max = max_index
        #print(l_max)
        #print(index)
        plt.figure(mof)
        plt.plot(p[:l_max],x[:l_max],label = ('Temp = %i K' % T))
        plt.xlim(0,100)
        plt.legend(loc = 'best')

        value = 5
        idx = (np.abs(p[:l_max]-value)).argmin()

        print(('Temp = %i K' % T), 'uptake @ 100 bar', '%0.3g' % np.max(x[:l_max]), 'mmol/g')
        print(('Temp = %i K' % T), 'uptake @ %0.3g bar' % p[idx], '%0.3g' % x[idx], 'mmol/g')

# matplotlib argument for showing plots
if 'noshow' not in sys.argv:
    plt.show()
