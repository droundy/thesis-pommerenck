#!/usr/bin/env python

from __future__ import division
import sys, os
import numpy as np
import readnew
from glob import glob
import matplotlib.pyplot as plt
import yaml
import os.path
import time # Need to wait some time if file is being written
import argparse

parser = argparse.ArgumentParser(description =
    """
        This script functions by reading in yaml files from a given directory
        and comparing with a given reference .DAT file. The heat capacity
        is computed and compared to the reference.
        
        EXAMPLE USAGE:
        python3 scripts-python-make-data/ising-multi-heat-capacity.py
        --file_dir=../ising-cp-data/
        --reference=../deft/papers/histogram/data/ising-32-reference-lndos.dat
        --save_dir=ising/data/comparison/N32
        --filename ising-sad-32
        --N=32
        --Emin=2048
        --Emax=0
        --seed_avg=8
    """)

parser.add_argument('--file_dir', required=True,
                    help='The directory where the files are located. \
                          Example:/home/jordan/ising-cp-data/')
parser.add_argument('--reference', required=True,
                    help='The directory where the reference file is located. \
                          Exmaple:data/ising-32-reference-lndos.dat')
parser.add_argument('--save_dir', required=True,
                    help='The directory where the data will be saved. \
                          Exmaple:data/comparison/N32')
parser.add_argument('--N', type=int, required=True,
                    help='The length of the ising model. Used to divide \
                          the number of moves.')
parser.add_argument('--Emin', type=int, required=True,
                    help='The energy at the entropy minimum.')
parser.add_argument('--Emax', type=int, required=True,
                    help='The energy at the entropy maximum.')
parser.add_argument('--seed_avg', type=int, required=True,
                    help='The number of seed files to average over.')
parser.add_argument('--filename', required=True, nargs='+',
                    help='The filename to perform analysis on. \
                          Exmaple:ising-sad-32')

args = parser.parse_args()

# Rename all argparse parameters.
file_dir = args.file_dir
ref = args.reference
save_dir = args.save_dir
N = args.N
Emin = args.Emin
Emax = args.Emax
seed_avg = args.seed_avg

# The correct way to handle accepting multiple arguments using nargs
# to create a list.
# '+' == 1 or more.
# '*' == 0 or more.
# '?' == 0 or 1.
filename = args.filename

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

T = np.arange(1.5, 5, 0.01)

# Compute the reference energy, entropy, and heat capacity once.
try:
    eref, lndosref, Nrt_ref = readnew.e_lndos_ps(ref)
except:
    eref, lndosref = readnew.e_lndos(ref)

#cvref = heat_capacity(T, eref, lndosref)
cvref = heat_capacity(T, eref[0:Emin-Emax+1], lndosref[0:Emin-Emax+1])

for f in filename:
    err_in_S = []
    err_max = []
    min_moves = []
    
    err_in_cv = []
    name = '%s.yaml' % (f)
    for n in range(1, seed_avg+1):
        #try:
            name = '%s-s%s.yaml' % (f, n)
            print(('trying filename ', name))

            # Read YAML file
            if os.path.isfile(file_dir + name):
                with open(file_dir + name, 'r') as stream:
                    yaml_data = yaml.load(stream)
            else:
                print(('unable to read file', file_dir + name))
                raise ValueError("%s isn't a file!" % (file_dir + name))

            data = yaml_data
            data['bins']['histogram'] = np.array(data['bins']['histogram'])
            data['bins']['lnw'] = np.array(data['bins']['lnw'])

            data['movies']['entropy'] = np.array(data['movies']['entropy'])
            lndos = data['movies']['entropy']
            energies = data['movies']['energy']
            N_save_times = len(data['movies']['entropy'])
            try:
                maxyaml = energies.index(-Emin)
            except:
                my_minE = energies[0]
                num_new_energies = int(my_minE - (-Emin))
                print(('num_new_maxyaml', num_new_energies))
                lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                lndos_new[:, num_new_energies:] = lndos[:,:]
                lndos = lndos_new
                energies = [0]*num_new_energies + energies
                maxyaml = 0

            try:
                minyaml = energies.index(-Emax)
            except:
                my_maxE = energies[-1]
                num_new_energies = -int(my_maxE - (-Emax))
                print(('num_new_minyaml', num_new_energies))
                lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                lndos_new[:, :lndos.shape[1]] = lndos[:,:]
                lndos = lndos_new
                minyaml = lndos.shape[1]-1

            errorinentropy = np.zeros(N_save_times)
            maxerror = np.zeros(N_save_times)
            cverror = np.zeros(N_save_times)

            for i in range(0, N_save_times):
                # below just set average S equal between lndos and lndosref
                if 'ising' in save_dir:

                    ising_norm = lndos[i][maxyaml:minyaml+1] # remove impossible state

                    ising_lndos = lndos[i][maxyaml:minyaml+1][::-1] # remove impossible state

                    # the states are counted backward hence the second to last state would be at index = 1
                    ising_norm = np.delete(ising_norm, [1])
                    ising_lndos = np.delete(ising_lndos, [len(ising_lndos)-2])

                    ising_E = np.array(energies[maxyaml:minyaml+1])
                    ising_E = np.delete(ising_E, [len(ising_E)-2])
                    
                    # invoke np.flip since ising_E is indexed backward!
                    # this is critical for my_cv_error or you get wrong answer.
                    flip_ising_E = np.flip(np.copy(ising_E))

                    norm_factor = np.mean(ising_norm) - np.mean(lndosref[0:Emin-Emax+1])
                    doserror = ising_lndos - lndosref[0:Emin-Emax+1] - norm_factor

                    my_cv_error = heat_capacity(T, flip_ising_E, ising_lndos) - cvref

                else:
                    print('Error! ising must be in save_dir pathname.')

                errorinentropy[i] = np.sum(abs(doserror))/len(doserror) #- np.mean(doserror)
                maxerror[i] = np.amax(doserror) - np.amin(doserror)
                cverror[i] = np.sum(abs(my_cv_error))/len(my_cv_error)
                max_cverror[i] = np.amax(abs(my_cv_error))
                # plt.figure('foo')
                # plt.clf()
                # plt.plot(T, heat_capacity(T, flip_ising_E, ising_lndos), label='myself')
                # plt.plot(T, cvref, label='cvref')
                # plt.legend(loc='best')
                # plt.title(f'save {i}')
                # plt.figure('foo diff')
                # plt.clf()
                # plt.plot(T, heat_capacity(T, flip_ising_E, ising_lndos) - cvref)
                # plt.title(f'error {i}')
                # plt.pause(0.000001)
                # DEBUG by looking into DOS and DOS ERROR!
                
            
            # remove N from moves in yaml file because N is added back in the
            # comparison-plot script
            moves = data['movies']['time']
            if len(min_moves) == 0 or len(min_moves) > len(moves):
                min_moves = np.array(moves)/N
            errorinentropy = errorinentropy[:len(moves)]
            maxerror = maxerror[:len(moves)]
            err_in_S.append(errorinentropy)
            err_max.append(maxerror)
            
            err_in_cv.append(max_cverror)


    for i in range(len(err_in_S)):
        err_in_S[i] = err_in_S[i][:len(min_moves)]
        err_max[i] = err_max[i][:len(min_moves)]
        
        err_in_cv[i] = err_in_cv[i][:len(min_moves)]

    errorinentropy = np.average(err_in_S, axis=0)
    maxmean = np.amax(err_in_S, axis=0)
    minmean = np.amin(err_in_S, axis=0)
    maxerror = np.average(err_max, axis=0)
    
    cverror = np.average(err_in_cv, axis = 0)
    max_mean_cv = np.amax(err_in_cv, axis=0)
    min_mean_cv = np.amin(err_in_cv, axis=0)

    dirname = '%s-%s' % (save_dir, name.replace('-s%s.yaml' %seed_avg, ''))
    
    print('saving to', dirname)
    try:
        os.mkdir(dirname)
    except OSError:
        pass
    else:
        print("Successfully created the directory %s " % dirname)
    np.savetxt('%s/errors.txt' %(dirname),
              np.c_[min_moves, errorinentropy, maxerror, minmean, maxmean, cverror, min_mean_cv, max_mean_cv],
              fmt = ('%.4g'),
              delimiter = '\t',
              header = 'iterations\t errorinentropy\t maxerror\t minmean\t maxmean\t cverror\t mincv\t maxcv\t(generated with python %s' % ' '.join(sys.argv))
