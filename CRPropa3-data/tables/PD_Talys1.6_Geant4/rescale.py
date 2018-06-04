from pylab import *
from xstools import xsGeant4

# -------------------------------------------------
# Rescale TALYS cross-sections to those of Kossov / Geant4
# -------------------------------------------------

# 1) Calculate the Geant4 cross-sections
eps = genfromtxt('eps.txt')  # photon energies [MeV]
idiv = eps.searchsorted(140)  # divide at eps = 140 MeV
eps_a = eps[:idiv]
eps_b = eps[idiv:]

# mean "extrapolation" for high energies from TALYS
# expo = (log10(mean(xs[:,499])) - log10(mean(xs[:,413]))) / (log10(eps[499]) - log10(eps[413]))
expo = -1.9

fsum = open('xs_sum.txt', 'w')
fmt = '%i\t%i' + '\t%g'*301 + '\n'

isotopes = genfromtxt('isotopes.txt', dtype=int)
xs = zeros((len(isotopes), 301))

for i,(Z,N,A) in enumerate(isotopes):
    xs_a = xsGeant4.crosssection(Z, N, eps_a)
    xs_b = xs_a[-1] * (eps_b / eps_a[-1])**expo
    xs[i] = r_[xs_a, xs_b]
    fsum.write(fmt % ((Z, N) + tuple(xs[i])))

fsum.close()


# 2) Rescale partial cross-sections from TALYS
d1 = genfromtxt('../PD_Talys1.6_Kahn/xs_sum.txt',
    dtype=[('Z',int), ('N',int), ('xs','301f8')])
d2 = genfromtxt('../PD_Talys1.6_Kahn/xs_thin.txt',
    dtype=[('Z',int), ('N',int), ('ch',int), ('xs','301f8')])

# scaling factor
r = xs / d1['xs']
r[isnan(r)] = 0
r[isinf(r)] = 0

for i,(Z,N,A) in enumerate(isotopes):
    idx = (d2['Z'] == Z) * (d2['N'] == N)
    d2['xs'][idx] *= r[i]

savetxt('xs_thin.txt', c_[d2['Z'], d2['N'], d2['ch'], d2['xs']],
    fmt='%i\t%i\t%i'+'\t%g'*301)
