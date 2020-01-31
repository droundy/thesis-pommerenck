#!/usr/bin/python3

import numpy as np
import os
from glob import glob
import argparse

parser = argparse.ArgumentParser(description =
    """
        This script functions by reading in energy and entropy data given
        the number 'N' forming the LJ cluster. The heat capacity is computed
        and compared to a benchmark (which is also read in).  A heat capacity
        error .TXT file is saved in the appropriate data directory.
    """)

parser.add_argument('--N', metavar='cluster_size', type=int, default=1,
                    help='The number of atoms in the LJ cluster -- Default = 1')
parser.add_argument('--data_dir', required=True,
                    help='The data directory to search through -- Required \
                          Example: /home/jordan/sad-monte-carlo/')
# parser.add_argument('--need_wide_cv', action='store_true',
#                     help='A boolean that determines whether we plot a narrow range --Default = False')

args = parser.parse_args()
# Possibly add in (3) more command line arguments that allow us to exactly
# specify the temperature range to plot.

N=args.N
#need_wide_cv = args.need_wide_cv
datadir = args.data_dir

def heat_capacity(T, E, S, Num):
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
    return C + Num*1.5

savedir = '../lj-cluster/data/lj%i/' % N
os.makedirs(savedir, exist_ok=True)

# Define the temperature range
T = np.arange(0.01, 0.05000001, 0.05000001*0.01)
wideT = np.arange(0.01, 0.4000001, 0.001)

CV = None
old_cvs = []

fnames = []

# Find all of the files with the matching number of particles N.
glob_names = glob((datadir + '*-%i-*' % N))
if glob_names == []:
    glob_names = glob((datadir + '*-%i.*' % N))

for glob in glob_names:
    fnames.append(glob.split('.')[0].replace(datadir,''))

fnames = list(set(fnames.copy()))

# Sort the benchmark files and put them first in the list for comparison.
bench_files = []
for i in range(len(fnames)):
    if 'bench' in fnames[i]:
        bench_files.append(fnames[i])

bench_files = sorted(bench_files)
fnames = bench_files + [x for x in fnames if x not in bench_files]

for fname in fnames:
    print(fname)
    time = np.loadtxt(datadir+fname+'.time')
    my_energy = np.loadtxt(datadir+fname+'.energy')
    my_entropy = np.loadtxt(datadir+fname+'.entropy')

    if CV is None:
        Ebest = my_energy;
        Sbest = my_entropy[-1,:]
        CV = heat_capacity(T, Ebest, Sbest, N)
        if 'bench' in fname:
            np.savetxt(savedir +'bench-cv.txt',
                    np.array([T, CV]).transpose(),
                    fmt='%.4g', # good enough for our plot
                    );

    cv_error = []
    cv_max_error = []
    myt = []

    for t in range(len(my_entropy[:, 0])):
        if time[t] < 1e3:
            continue
        myt.append(time[t])
        mycv = heat_capacity(T, my_energy, my_entropy[t,:], N)

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
    if 'bench' not in fname:
        np.savetxt(savedir+os.path.basename(fname)+'-cv-error.txt',
                np.array([myt, cv_error, cv_max_error]).transpose(),
                fmt='%.3g', # low resolution is good enough for error data.
                )

import matplotlib.pyplot as plt

plt.figure()
plt.loglog(myt,cv_error)

plt.figure()
plt.plot(wideT,heat_capacity(wideT, Ebest, Sbest, N))
plt.show()
