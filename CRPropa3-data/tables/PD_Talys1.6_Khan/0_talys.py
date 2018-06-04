from numpy import *
import os
import time


# list of isotopes (Z <= 26, A <= 30, lifetime > )
isotopes = genfromtxt('isotopes.txt')
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


# steering card
s  = 'projectile g  \n' # incident photon
s += 'energy eps.txt\n' # file with incident photon energies
s += 'element %s    \n' # target element
s += 'mass %i       \n' # target mass number
s += 'strength 1    \n' # GDR parameterization (E1 strength function): 1 = Kopecky-Uhl / generalized Lorentzian
s += 'gnorm 1.0     \n' # GDR normalization factor
s += 'channels y    \n' # write exclusive cross-sections output files
s += 'maxchannel 8  \n' # maximum number of emitted particles per for which output files are written
s += 'isomer 1.e38  \n' # don't consider isomers separately
s += 'fileresidual n\n' # don't write residual production output files


# incident photon energies: 0.2 - 200 MeV in steps of dlogE = 0.01
eps = logspace(log10(0.2), log10(200), 301)
savetxt('eps.txt', eps, fmt='%.6g')



for (Z,N,A) in isotopes[:2]:
    print Z, N
    time.sleep(3)

    folder = '%i-%i' % (Z, A)
    try:
        os.mkdir(folder)
    except:
        pass

    os.chdir(folder)

    f = open('input', 'w')
    f.write(input_str % (elements[Z], A))
    f.close()

    os.system('talys < input > output')
    os.chdir('..')
