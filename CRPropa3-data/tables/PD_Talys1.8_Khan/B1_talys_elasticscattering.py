from numpy import *
import os
import time


elements = {
    1  :  'H',
    2  : 'He',
    3  : 'Li',
    4  : 'Be',
    5  :  'B',
    6  :  'C',
    7  :  'N',
    8  :  'O',
    9  :  'F',
    10 : 'Ne',
    11 : 'Na',
    12 : 'Mg',
    13 : 'Al',
    14 : 'Si',
    15 :  'P',
    16 :  'S',
    17 : 'Cl',
    18 : 'Ar',
    19 :  'K',
    20 : 'Ca',
    21 : 'Sc',
    22 : 'Ti',
    23 :  'V',
    24 : 'Cr',
    25 : 'Mn',
    26 : 'Fe'}


# Steering card
s  = 'projectile g  \n' # incident photon
s += 'energy %s     \n' # file with incident photon energies
s += 'element %s    \n' # target element
s += 'mass %i       \n' # target mass number
s += 'strength 1    \n' # GDR parameterization (E1 strength function): 1 = Kopecky-Uhl / generalized Lorentzian
s += 'gnorm 1.0     \n' # GDR normalization factor
s += 'channels n    \n' # don't write exclusive cross-sections output files
s += 'maxchannel 8  \n' # maximum number of emitted particles per for which output files are written
s += 'isomer 1.e38  \n' # don't consider isomers separately
s += 'fileresidual n\n' # don't write residual production output files


# Create file of incident photon energies: 0.002 - 200 MeV in steps of dlogE = 0.01
eps = logspace(log10(0.002), log10(200), 501)
savetxt('eps_elastic.txt', eps, fmt='%.6g')
fname_eps = os.getcwd() + '/eps_elastic.txt'


# Run TALYS for each isotope (Z <= 26, A <= 30, lifetime > 2s)
# Note: TALYS cannot simulate H, He as targets and fails for Li and Be
for (Z,N,A) in genfromtxt('isotopes.txt'):
    print Z, N
    time.sleep(3)

    folder = '%i-%i/elastic_scattering' % (Z, A)
    try:
        os.mkdir(folder)
    except:
        pass

    os.chdir(folder)

    f = open('input', 'w')
    f.write(s % (fname_eps, elements[Z], A))
    f.close()

    os.system('talys < input > output')
    os.chdir('../..')
