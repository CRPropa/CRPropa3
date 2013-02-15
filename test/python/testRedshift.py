#!/usr/bin/env python
from mpc import *
from pylab import *


# Redshift - comoving distance relation
redshift = Redshift()

N = 100
Z = logspace(-4, 2, N)
D1 = zeros(N)
for i in range(N):
	D1[i] = redshift.redshift2ComovingDistance(Z[i]) / Mpc

H0 = redshift.hubbleRate(0)
D2 = c_light * Z / H0 / Mpc

figure()
plot(Z, D1, label='Numerical Integration')
plot(Z, D2, 'r--', label='Small Redshift Approximation')
legend(loc='upper left', frameon=0)
xlabel('Redshift z')
ylabel('Comoving Distance [Mpc]')
loglog()
grid()
savefig('Redshift.png', bbox_inches='tight')
show()

