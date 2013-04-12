#!/usr/bin/env python
from mpc import *
from pylab import *


setCosmologyParameters(0.71, 0.27)

N = 100
Z = logspace(-4, 2, N)
Dc = zeros(N) # comoving distance
Dl = zeros(N) # comoving distance
Dt = zeros(N) # light travel distance
for i in range(N):
    Dc[i] = redshift2ComovingDistance(Z[i]) / Mpc
    Dl[i] = redshift2LuminosityDistance(Z[i]) / Mpc
    Dt[i] = redshift2LightTravelDistance(Z[i]) / Mpc

H0 = hubbleRate(0)
Dlin = c_light * Z / H0 / Mpc # Hubble's law

figure()
plot(Z, Dc, 'b', label="Comoving distance")
plot(Z, Dl, 'g', label="Luminosity distance")
plot(Z, Dt, 'r', label="Light travel distance")
plot(Z, Dlin, 'k--', label="Naive Hubble")

legend(loc='upper left', frameon=0)
xlabel('Redshift z')
ylabel('Distance [Mpc]')
loglog()
grid()
savefig('Redshift.png', bbox_inches='tight')
show()

