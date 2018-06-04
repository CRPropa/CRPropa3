from numpy import *
import os

# -----------------------------------------------
# Collect all exclusive cross sections from TALYS
# -----------------------------------------------

def isBogus(Z, N, channel):
    s = channel  # expect a string
    dZ = int(s[1]) + int(s[2]) + int(s[3]) + 2*int(s[4]) + 2*int(s[5])
    dN = int(s[0]) + int(s[2]) + int(s[3]) +   int(s[4]) + 2*int(s[5])
    if (dZ + dN) == 0:
        print '    no photo-disintegration:', channel
        return True
    if (dZ + dN) == (Z + N):
        print '    no nucleon left:', channel
        return True
    if dZ > dZ:
        print '    too many protons lost:', channel
        return True
    if dN > dN:
        print '    too many neutrons lost:', channel
        return True
    return False


fexcl = open('xs_excl.txt', 'w')
fexcl.write('# Z\tN\tchannel\txs\n')
fexcl.write('#cross sections [mb] for incident photon energies eps = 0.2 - 200 MeV in steps of logE = 0.01\n')
fmt = '%i\t%i\t%s' + '\t%.4g'*301 + '\n'  # output format

isotopes = genfromtxt('isotopes.txt', dtype=int)

for (Z,N,A) in isotopes:
    print Z, N
    folder = '%i-%i/' % (Z, A)

    for f in os.listdir(folder):
        if not( f.startswith('xs') and f.endswith('tot') ):
            continue  # consider the xs....tot files only

        channel = f.strip('xs.tot')
        if isBogus(Z, N, channel):
            continue  # skip bogus channels

        xs = genfromtxt(folder + f, usecols=1)
        fexcl.write( fmt % ((Z, N, channel) + tuple(xs)) )

fexcl.close()
