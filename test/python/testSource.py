from mpc import *
from pylab import *

source = CompositeSource(Vector3d(0.), 10, 100, -1)

def add(A, Z, w):
    source.addToComposition(getNucleusId(A, Z), w * A)

# DuVernois
add(1, 1, 92000)
add(4, 2, 13000)
add(12, 6, 447.4)
add(14, 7, 34.2)
add(16, 8, 526.3)
add(20, 10, 58.0)
add(23, 11, 3.2)
add(24, 12, 108.0)
add(27, 13, 7.8)
add(28, 14, 100)
add(32, 16, 13.1)
add(40, 18, 2.2)
add(40, 20, 6.5)
add(52, 24, 1.5)
add(55, 25, 1.1)
add(56, 26, 97)

state = ParticleState()
N = 10000
Etot = zeros(N)
E1, E2, E3, E4 = [], [], [], []
for i in range(N):
    source.prepare(state)
    lE = log10(state.getEnergy()) + 18
    Z = state.getChargeNumber()
    Etot[i] = lE
    if Z == 1:
        E1.append(lE)
    elif Z == 2:
        E2.append(lE)
    elif Z <= 8:
        E3.append(lE)
    else:
        E4.append(lE)

E1 = array(E1)
E2 = array(E2)
E3 = array(E3)
E4 = array(E4)

figure()
hist(E1, weights=ones(len(E1))*1./N, bins=20, range=[19, 21.5], histtype='step', lw=2, label='H')
hist(E2, weights=ones(len(E2))*1./N, bins=20, range=[19, 21.5], histtype='step', lw=2, label='He')
hist(E3, weights=ones(len(E3))*1./N, bins=20, range=[19, 21.5], histtype='step', lw=2, label='$Z \leq 8$')
hist(E4, weights=ones(len(E4))*1./N, bins=20, range=[19, 21.5], histtype='step', lw=2, label='$Z > 8$')
ylim(1e-4, 1)
xlabel('$\log_{10}(E/\mathrm{eV})$')
ylabel('Probability')
semilogy()
grid()
legend()
savefig('Source.png')
