#!/usr/bin/env python

from __future__ import division
import sys, os, matplotlib
import numpy as np
import readnew
from glob import glob
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')
matplotlib.rcParams['figure.figsize'] = (5, 4)

if 'noshow' in sys.argv:
        matplotlib.use('Agg')

import yaml
import os.path
import time # Need to wait some time if file is being written
import argparse
import colors
import pandas as pd

parser = argparse.ArgumentParser(description =
    """
        This script functions by reading in a .CSV file using
        Pandas and plots the data.
        
        EXAMPLE USAGE:
        python3 scripts-python/ising-cv-plot.py
        --file_dir=../ising-cp-data/
        --save_dir=ising/data/heat-capacity
        --filename ising-sad-32
        --N=32
    """)

parser.add_argument('--file_dir', required=True,
                    help='The directory where the files are located. \
                          Example:/home/jordan/ising-cp-data/')
# parser.add_argument('--save_dir', required=True,
#                     help='The directory where the data will be saved. \
#                           Exmaple:data/comparison/N32')
parser.add_argument('--N', type=int, required=True,
                    help='The length of the ising model. Used to divide \
                          the number of moves.')
# parser.add_argument('--filename', required=True,
#                     help='The filename to perform analysis on. \
#                           Exmaple:ising-sad-32')

args = parser.parse_args()

# Rename all argparse parameters.
file_dir = args.file_dir
# save_dir = args.save_dir
N = args.N
# filename = args.filename

cv_data = pd.read_csv('%s/N%s-heat-capacity.csv' % (file_dir,N),delimiter='\t',encoding='utf-8',engine='python')

#cv_data['cvref'] = cv_data['cvref'].astype(float)
print(cv_data.head(10))

cv_headers = list(cv_data)[1:]
print(cv_headers)

# Begin plotting the heat capacity
plt.figure('heat capacity plot')
Temp = np.array(pd.to_numeric(cv_data['Temperature'], errors='coerce'))
for name in cv_headers:
    cv = pd.to_numeric(cv_data[name], errors='coerce')
    colors.plot(1/Temp, cv / N**2, method=name)

colors.legend(loc='best')
plt.xlabel(r'$k_B\beta$')
plt.ylabel(r'$c_V$ / $ k_B$')
plt.xlim(0.4,0.5)
if N == 32:
    plt.ylim(0.6,2.4)

plt.savefig('../ising/N%i-Cv.pdf' % N)

plt.show()