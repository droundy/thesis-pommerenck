from __future__ import division, print_function
import sys, os, matplotlib
import numpy as np
from collections import OrderedDict

matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')
matplotlib.rcParams['figure.figsize'] = (5, 4)

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
import colors

if os.path.exists('ising/data'):
    os.chdir('..')
    
if os.path.exists('square-well/data'):
    os.chdir('..')

energy = int(sys.argv[1])
filebase = sys.argv[2]
transcale = sys.argv[3]
system = sys.argv[4]

tex_filebase = filebase.replace('.','_') # latex objects to extra "." characters

methods = ['-sad3', '-sad3-s1', '-sad3-s2',
            '-tmmc', '-tmi', '-tmi2', '-tmi3', '-toe', '-toe2', '-toe3',
            '-vanilla_wang_landau']
if 'allmethods' not in sys.argv:
    methods = ['-sad3','-wl-50','-wl-inv-t-50','-wl-inv-t-256','-sad-50',
               '-sad-256','-wl-256','-ising-sad-32','-ising-wl-32','-ising-wl-inv-t-32',
               '-ising-sad-128','-ising-wl-128','-ising-wl-inv-t-128',
               '-ising-wl-32-minGamma', '-ising-wl-128-minGamma', '-ising-wl-128-minGamma-1e6']#, '-ising-sad-128-T15'] #This file has a lot of data
    # methods = [m+s for m in methods for s in '', '-s1', '-s2', '-s3', '-s4']
    if transcale == 'slow':
        methods = ['-sad3-slow','-sad-slow-T13','-wl-50-slow','-wl-inv-t-50-slow','-sad-50-slow',
                   '-tmmc-slow', '-vanilla_wang_landau-slow','-vanilla_wang_landau-T13-slow','-sad3-T13-slow']
    if transcale == 'fast':
        methods = ['-sad3-fast','-tmmc-fast', '-vanilla_wang_landau-fast']

# For WLTMMC compatibility with LVMC
lvextra = glob('data/comparison/%s-wltmmc*' % filebase)
split1 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra]
split2 = [i.split('-m', 1)[0] for i in split1]
for meth in split2:
    if "default" in transcale:
        if meth[-3:] != '-tm' and "slow" not in meth and "fast" not in meth:
            methods.append('-%s' % meth)
    if "slow" in transcale:
        if meth[-3:] != '-tm' and "slow" in meth:
            methods.append('-%s' % meth)
    if "fast" in transcale:
        if meth[-3:] != '-tm' and "fast" in meth:
            methods.append('-%s' % meth)

# For SAMC compatibility with LVMC
if system == 's000':
    lvextra1 = glob('data/comparison/%s/%s-samc*' % (system,filebase))
else:
    lvextra1 = glob('data/comparison/%s-*samc*' % filebase)
split3 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra1]
split4 = [i for i in split3 if i[-3:-1] != '-s']

for meth in split4:
    if "default" in transcale:
        if meth[-3:] != '-tm' and "slow" not in meth and "fast" not in meth:
            methods.append('-%s' % meth)
    if "slow" in transcale:
        if meth[-3:] != '-tm' and "slow" in meth:
            methods.append('-%s' % meth)
    if "fast" in transcale:
        if meth[-3:] != '-tm' and "fast" in meth:
            methods.append('-%s' % meth)


best_ever_max = 1e100
max_time = 0
min_error = 1e200
print('methods are', methods)
for method in methods:
    print('trying method', method)
    try:
        if system == 's000':
                dirname = 'data/comparison/%s/%s%s/' % (system,filebase,method)
        else:
                dirname = 'data/comparison/%s%s/' % (filebase,method)

        if not os.path.exists(dirname) or os.listdir(dirname) == []:
                continue

        if energy > 0:
                Nrt_at_energy, erroratenergy = np.loadtxt(dirname + 'energy-%s.txt' % energy, delimiter = '\t', unpack = True)
        data = np.loadtxt(dirname + 'errors.txt', delimiter = '\t', unpack = True)
        iterations = data[0]
        errorinentropy = data[1]
        maxerror = data[2]
        best_ever_max = min(best_ever_max, maxerror.min())

        #try:
            #my_e, my_S_error_vs_e = np.loadtxt(dirname + 'error-vs-energy.txt', delimiter = '\t', unpack = True)
            #plt.figure('error-vs-energy')
            #colors.plot(my_e, my_S_error_vs_e, method=method[1:])
        #except:
            #pass

        howManyPoints = 300 # We are definately using C++ code!
        if len(iterations) > howManyPoints:
                temp_iter = []
                temp_ee = []
                temp_me = []
                last_power = 0
                movie_time = 10.0**(1.0/8)
                for i in range(len(iterations)):
                        if iterations[i] >= movie_time**last_power:
                                temp_iter.append(iterations[i])
                                temp_ee.append(errorinentropy[i])
                                temp_me.append(maxerror[i])
                                while iterations[i] >= movie_time**last_power:
                                        last_power += 1
                iterations = np.array(temp_iter)
                errorinentropy = np.array(temp_ee)
                maxerror = np.array(temp_me)

        if system == 's000':
                N = filebase.split('-N')[-1]
                ff = filebase.split('-ff')[-1].split('-N')[0]
                ff = float(ff)
                # Get N directly from title.
                moves = iterations * float(N)
        if system == 'ising':
                N = filebase.split('N')[-1]
                moves = iterations * float(N)

        max_time = max(max_time, moves.max())

        if energy > 0:
                plt.figure('error-at-energy-iterations')
                colors.plot(moves, erroratenergy, method=method[1:])
                plt.title('Error at energy %g %s' % (energy,filebase))
                plt.xlabel('# moves')
                plt.ylabel('error')
                colors.legend()
                plt.savefig('figs/%s-error-energy-%g.pdf' % (tex_filebase, energy), transparent=True)

                plt.figure('round-trips-at-energy' )
                colors.plot(moves, Nrt_at_energy, method = method[1:])
                plt.title('Round Trips at energy %g, %s' % (energy,filebase))
                plt.xlabel('# moves')
                plt.ylabel('Round Trips')
                colors.legend()
                plt.savefig('figs/%s-round-trips-%g.pdf' % (tex_filebase, energy), transparent=True)

                plt.figure('error-at-energy-round-trips')
                colors.plot(Nrt_at_energy[Nrt_at_energy > 0], erroratenergy[Nrt_at_energy > 0],
                         method = method[1:])
                plt.title('Error at energy %g %s' % (energy,filebase))
                plt.xlabel('Round Trips')
                plt.ylabel('Error')
                colors.legend()
                plt.savefig('figs/%s-error-energy-Nrt-%g.pdf' % (tex_filebase, energy), transparent=True)

        plt.figure('maxerror')
        colors.loglog(moves, maxerror, method = method[1:])
        plt.xlabel(r'Moves')
        plt.ylabel(r'Maximum Entropy Error')
        #plt.title('Maximum Entropy Error vs Iterations, %s' %filebase)
        colors.legend()

        if type(iterations) is not np.float64:
                plt.figure('errorinentropy')
                if data.shape[0] > 4:
                    minmean = data[3]
                    maxmean = data[4]
                    plt.fill_between(moves, minmean, maxmean,
                                     edgecolor='none', linewidth=0,
                                        color=colors.color(method[1:]),
                                        alpha=0.1, zorder=-51)
                my_S_error = errorinentropy[0:len(iterations)]
                min_error = min(min_error, my_S_error[my_S_error > 0].min())
                colors.loglog(moves, my_S_error, method = method[1:])
                plt.xlabel(r'$\textrm{Moves}$')
                plt.ylabel(r'$\textrm{Average Entropy Error}$')
                #plt.title('Average Entropy Error at Each Iteration, %s' %filebase)
                if filebase.startswith('s000'):
                    if "default" in transcale:
                        plt.title(r'$N=%d$, $\eta = %g$ $\delta_0 = 0.05$' % (int(N), ff))
                    elif "slow" in transcale:
                        plt.title(r'$N=%d$, $\eta = %g$ $\delta_0 = 0.005$' % (int(N), ff))
                    elif "fast" in transcale:
                        plt.title(r'$N=%d$, $\eta = %g$ $\delta_0 = 0.5$' % (int(N), ff))
                else: # We are doing Ising not square well
                    if "default" in transcale:
                        plt.title(r'$N=%d$, $\delta_0 = 0.05$' % (int(N)))
                    elif "slow" in transcale:
                        plt.title(r'$N=%d$, $\delta_0 = 0.005$' % (int(N)))
                    elif "fast" in transcale:
                        plt.title(r'$N=%d$, $\delta_0 = 0.5$' % (int(N)))
                colors.legend()
                if data.shape[0] > 5:
                    plt.figure('heat_capacity')
                    if filebase == 'N32':
                        plt.xlim(1e5,1e12)
                        plt.ylim(1e-3,1e1)
                        cv_factor = (32*32)
                    else:
                        plt.xlim(1e8,10**(13.8))
                        plt.ylim(1e-3,1e1)
                        cv_factor = 128*128
                    try:
                        heat_capacity = data[5] / cv_factor
                        minmean_cv = data[6] / cv_factor
                        maxmean_cv = data[7] / cv_factor
                        plt.fill_between(moves, minmean_cv, maxmean_cv,
                                        edgecolor='none', linewidth=0,
                                        color=colors.color(method[1:]),
                                        alpha=0.1, zorder=-51)
                    except:
                        pass
                    colors.loglog(moves, heat_capacity, method = method[1:])
                    for i in np.arange(-8, 19, 1.0):
                        colors.loglog(moves, 10**i/np.sqrt(0.1*moves), method = r'1/sqrt(t)')
                    plt.xlabel(r'$\textrm{Moves}$')
                    plt.ylabel(r'Maximum error in $c_V$')
                    if filebase == 'N128':
                      colors.legend(loc='lower left')
                    else:
                      colors.legend()
                    plt.savefig('%s-Cv-error-%s.pdf' % (tex_filebase,transcale), transparent=True)

    except:
        raise
plt.figure('maxerror')
colors.loglog([1e6,max_time], [best_ever_max*np.sqrt(max_time/1e6), best_ever_max], method = r'$\frac{1}{\sqrt{t}}$')
colors.legend()
#plt.savefig('%s-max-entropy-error-%s.pdf' % (tex_filebase,transcale))

plt.figure('errorinentropy')
if filebase == 'N128':
    moves = np.array([1e5, 10**(13.8)])
else:
    moves = np.array([1e5, 1e12])
#colors.loglog(moves, min_error*np.sqrt(moves.max())/np.sqrt(moves), method = r'1/sqrt(t)')
for i in np.arange(-8, 19, 1.0):
    colors.loglog(moves, 10**i/np.sqrt(0.1*moves), method = r'1/sqrt(t)')
plt.xlim(moves[0], moves[1])
if filebase == 'periodic-ww1.30-ff0.30-N50':
    plt.ylim(1e-3, 1e4)
    if "slow" in transcale:
        plt.ylim(1e-2, 1e5)
elif filebase == 'periodic-ww1.30-ff0.30-N500':
    plt.ylim(1e-1, 1e3)
elif filebase == 'periodic-ww1.50-ff0.17-N256':
    plt.ylim(1e-3, 1e3)
elif filebase == 'N32':
    plt.ylim(1e-3,1e3)
elif filebase == 'N128':
    plt.ylim(1e-3,1e4)


colors.legend()
plt.tight_layout()
print('filename', '%s-entropy-error-%s.pdf' % (tex_filebase,transcale))
plt.savefig('%s-entropy-error-%s.pdf' % (tex_filebase,transcale))

#plt.figure('error-vs-energy')
#plt.ylim(-1,1)
#colors.legend()
#plt.tight_layout()
#plt.savefig('figs/%s-error-vs-energy-%s.pdf' % (tex_filebase,transcale))

plt.show()
