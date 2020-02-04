#!/usr/bin/python3

import numpy as np
import os
import yaml
import argparse

parser = argparse.ArgumentParser(description =
    """
        This script functions by reading in energy and entropy data given
        the number 'N' forming the N*N Ising model. The heat capacity is computed
        and compared to a benchmark (which is also read in).  A heat capacity
        error .TXT file is saved in the appropriate data directory.
    """)

parser.add_argument('--N', metavar='ising_length', type=int, default=1,
                    help='The length of the ising model -- Default = 1')
parser.add_argument('--data_dir', required=True,
                    help='The data directory to search through -- Required \
                          Example: /home/jordan/sad-monte-carlo/')
parser.add_argument('--fname', required=True,
                    help='The file name -- Required \
                          Example: ising-sad-32-s1.yaml')
parser.add_argument('--refname', required=True,
                    help='The reference file name -- Required \
                          Example: ising-sad-32-s1.yaml')
# parser.add_argument('--need_wide_cv', action='store_true',
#                     help='A boolean that determines whether we plot a narrow range --Default = False')

args = parser.parse_args()

N=args.N
#need_wide_cv = args.need_wide_cv
datadir = args.data_dir

def heat_capacity(T, E, S):
    """
        This function calculates the heat capacity by taking in
        the temperature 'T', energy 'E', and entropy 'S' as NumPy
        arrays.  The number 'Num' is also necessary to complete
        the heat capacity calculation.
    """
    C = np.zeros_like(T)
    for i in range(len(T)):
        boltz_arg = S - E/T[i]
        P = np.exp(boltz_arg - boltz_arg.max())
        P = P/P.sum()
        U = (E*P).sum()
        C[i] = ((E-U)**2*P).sum()/T[i]**2
    return C

savedir = '../ising/data/ising-%i/' % N
os.makedirs(savedir, exist_ok=True)

# Define the temperature range
T = np.arange(0.01, 0.05000001, 0.05000001*0.01)
wideT = np.arange(0.01, 5, 0.001)

old_cvs = []

print(('trying filename ', args.fname))

# Read YAML file
if os.path.isfile(args.data_dir + args.fname):
    with open(args.data_dir + args.fname, 'r') as stream:
        data = yaml.load(stream)
else:
    print(('unable to read file', args.data_dir + args.fname))
    raise ValueError("%s isn't a file!" % (args.data_dir + args.fname))

# Read YAML Reference file
if os.path.isfile(args.data_dir + args.refname):
    with open(args.data_dir + args.refname, 'r') as stream:
        ref_data = yaml.load(stream)
else:
    print(('unable to read file', args.data_dir + args.refname))
    raise ValueError("%s isn't a file!" % (args.data_dir + args.refname))

time = data['movies']['time']
my_energy  = data['movies']['energy']
my_entropy = data['movies']['entropy']

ref_energy  = ref_data['movies']['energy']
ref_entropy = ref_data['movies']['entropy'][-1]

Ebest = ref_energy
Sbest = ref_entropy
CV = heat_capacity(T, Ebest, Sbest)

cv_error = []
cv_max_error = []
myt = []

for t in range(len(my_entropy)):
    # if time[t] < 1e3:
        # continue
    myt.append(time[t])
    mycv = heat_capacity(T, my_energy, my_entropy[t])

    err = 0
    norm = 0
    for j in range(1, len(mycv)):
        err += abs(CV[j]-mycv[j])
        norm += 1.0
    cv_error.append(err/norm)
    cv_max_error.append(abs(CV-mycv).max())
    # if time[t] == 1e12:
    #     np.savetxt(datadir+os.path.basename(fname)+'-cv.txt',
    #                np.array([T, mycv]).transpose(),
    #                fmt='%.4g', # good enough for our plot
    #     );
    #     np.savetxt(datadir+os.path.basename(fname)+'-wide-cv.txt',
    #                np.array([wideT, heat_capacity(wideT, my_energy, my_entropy[t,:])]).transpose(),
    #                fmt='%.4g', # good enough for our plot
    #     );

#old_cvs.append(mycv)
# if 'bench' not in fname:
    # np.savetxt(savedir+os.path.basename(fname)+'-cv-error.txt',
            # np.array([myt, cv_error, cv_max_error]).transpose(),
            # fmt='%.3g', # low resolution is good enough for error data.
            # )

import matplotlib.pyplot as plt

plt.figure('cv-error')
plt.loglog(myt,cv_error)

plt.figure('heat-capacity')
plt.xlim(0.3,1)
plt.plot(np.reciprocal(wideT),heat_capacity(wideT, Ebest, Sbest))
plt.plot(np.reciprocal(wideT),heat_capacity(wideT, my_energy, my_entropy[-1]))
plt.show()
