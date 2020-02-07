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
parser.add_argument('--N', type=int,
                    help='The length of the ising model. Used to divide \
                          the number of moves.')
parser.add_argument('--Emin', type=int,
                    help='The energy at the entropy minimum.')
parser.add_argument('--Emax', type=int,
                    help='The energy at the entropy maximum.')


args = parser.parse_args()

file_dir = args.file_dir                # change name from filename_location
ref = args.reference                    # reference to ref
save_dir = args.save_dir                # filebase



# seed_avg = int(sys.argv[8])

# filename = sys.argv[9:]
# print(('filenames are ', filename))

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

#T = np.arange(0.01, 0.05000001, 0.05000001*0.01)
T = np.arange(0.01, 5, 0.001)


#n_energies = int(Emin - Emax+1)

try:
    eref, lndosref, Nrt_ref = readnew.e_lndos_ps(ref)
except:
    eref, lndosref = readnew.e_lndos(ref)

CV_ref = heat_capacity(T, eref, lndosref)

plt.plot(T,CV_ref)
plt.show()

# for f in filename:
    # err_in_S = []
    # err_max = []
    # min_moves = []
    # name = '%s.yaml' % (f)
    # for n in range(1, seed_avg+1):
        # #try:
            # name = '%s-s%s.yaml' % (f, n)
            # print(('trying filename ', name))

            # # Read YAML file
            # if os.path.isfile(filename_location + name):
                # with open(filename_location + name, 'r') as stream:
                    # yaml_data = yaml.load(stream)
            # else:
                # print(('unable to read file', filename_location + name))
                # raise ValueError("%s isn't a file!" % (filename_location + name))

            # data = yaml_data
            # data['bins']['histogram'] = np.array(data['bins']['histogram'])
            # data['bins']['lnw'] = np.array(data['bins']['lnw'])

            # data['movies']['entropy'] = np.array(data['movies']['entropy'])
            # lndos = data['movies']['entropy']
            # energies = data['movies']['energy']
            # N_save_times = len(data['movies']['entropy'])
            # try:
                # maxyaml = energies.index(-Emin)
            # except:
                # my_minE = energies[0]
                # num_new_energies = int(my_minE - (-Emin))
                # print(('num_new_maxyaml', num_new_energies))
                # lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                # lndos_new[:, num_new_energies:] = lndos[:,:]
                # lndos = lndos_new
                # energies = [0]*num_new_energies + energies
                # maxyaml = 0

            # try:
                # minyaml = energies.index(-Emax)
            # except:
                # my_maxE = energies[-1]
                # num_new_energies = -int(my_maxE - (-Emax))
                # print(('num_new_minyaml', num_new_energies))
                # lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                # lndos_new[:, :lndos.shape[1]] = lndos[:,:]
                # lndos = lndos_new
                # minyaml = lndos.shape[1]-1

            # try:
                # eref, lndosref, Nrt_ref = readnew.e_lndos_ps(ref)
            # except:
                # eref, lndosref = readnew.e_lndos(ref)

            # errorinentropy = np.zeros(N_save_times)
            # maxerror = np.zeros(N_save_times)

            # for i in range(0, N_save_times):
                # # below just set average S equal between lndos and lndosref
                # if TRUE:
                    # # if using yaml as a reference the range is from 0 to len while for C++ the range is
                    # # from maxref to minref + 1
                    # if 'ising' in filebase:

                        # ising_norm = lndos[i][maxyaml:minyaml+1] # remove impossible state

                        # ising_lndos = lndos[i][maxyaml:minyaml+1][::-1] # remove impossible state

                        # # the states are counted backward hence the second to last state would be at index = 1
                        # ising_norm = np.delete(ising_norm, [1])
                        # ising_lndos = np.delete(ising_lndos, [len(ising_lndos)-2])

                        # ising_E = np.array(energies[maxyaml:minyaml+1])
                        # ising_E = np.delete(ising_E, [len(ising_E)-2])


                        # norm_factor = np.mean(ising_norm) - np.mean(lndosref[0:minref-maxref+1])
                        # doserror = ising_lndos - lndosref[0:minref-maxref+1] - norm_factor
                    # else:
                        # norm_factor = np.mean(lndos[i][maxyaml:minyaml+1]) - np.mean(lndosref[0:minref-maxref+1])
                        # doserror = lndos[i][maxyaml:minyaml+1][::-1] - lndosref[0:minref-maxref+1] - norm_factor
                # else:
                    # norm_factor = np.mean(lndos[i][maxyaml:minyaml+1]) - np.mean(lndosref[maxref:minref+1])
                    # doserror = lndos[i][maxyaml:minyaml+1][::-1] - lndosref[maxref:minref+1] - norm_factor

                # errorinentropy[i] = np.sum(abs(doserror))/len(doserror) #- np.mean(doserror)
                # maxerror[i] = np.amax(doserror) - np.amin(doserror)

            # # remove N from moves in yaml file because N is added back in the
            # # comparison-plot script
            # moves = data['movies']['time']
            # if min_moves == [] or len(min_moves) > len(moves):
                # min_moves = np.array(moves)/N
            # errorinentropy = errorinentropy[:len(moves)]
            # maxerror = maxerror[:len(moves)]
            # err_in_S.append(errorinentropy)
            # err_max.append(maxerror)


    # for i in range(len(err_in_S)):
        # err_in_S[i] = err_in_S[i][:len(min_moves)]
        # err_max[i] = err_max[i][:len(min_moves)]
    # errorinentropy = np.average(err_in_S, axis=0)
    # maxmean = np.amax(err_in_S, axis=0)
    # minmean = np.amin(err_in_S, axis=0)
    # maxerror = np.average(err_max, axis=0)
    # dirname = 'data/comparison/%s-%s' % (filebase, name.replace('-s%s.yaml' %seed_avg, ''))
