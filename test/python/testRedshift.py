#!/usr/bin/env python
from mpc import *
import numpy as np
import matplotlib.pyplot as plt



redshift = Redshift()

# Redshift - comoving distance relation
N = 100
Z = np.logspace(-4, 2, N)
D1 = np.zeros(N)
for i in range(N):
	D1[i] = redshift.getDistance(Z[i]) / Mpc

H0 = redshift.getHubbleRate(0)
D2 = c_light * Z / H0 / Mpc

plt.figure()
plt.plot(Z, D1, label='Numerical Integration')
plt.plot(Z, D2, 'r--', label='Small Redshift Approximation')
plt.legend(loc='upper left', frameon=0)
plt.xlabel('Redshift z')
plt.ylabel('Comoving Distance [Mpc]')
plt.loglog()
plt.grid()
plt.savefig('Redshift.png', bbox_inches='tight')
plt.show()

