from numpy import *
import os


isotopes = genfromtxt('isotopes.txt')
path = './'  # path to simulated TALYS output

# -----------------------------------------------------------------------------
# Collect all cross sections for elastic scattering from TALYS
# -----------------------------------------------------------------------------
fexcl = open('xs_elastic.txt', 'w')
fexcl.write('# Z\tN\txs\n')
fexcl.write('#cross sections [mb] for incident photon energies eps = 0.002 - 200 MeV in steps of logE = 0.01\n')
fmt = '%i\t%i' + '\t%.4g'*501 + '\n'  # output format

for (Z,N,A) in isotopes:
    print Z, N
    folder = path + '%i-%i/elastic_scattering/' % (Z, A)
    if not (os.path.exists(folder + "elastic.tot")):
        print "no elastic scattering data"
        continue
    xs = genfromtxt(folder + "elastic.tot", usecols=1)
    fexcl.write( fmt % ((Z, N) + tuple(xs)) )

fexcl.close()
