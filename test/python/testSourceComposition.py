from mpc import *
from pylab import *

# Simulates the integrated relative abundance from a source that accelerates
# particles up to a maximum rigidity.
# min. energy = 10 EeV
# max. rigidity = 100 EeV / Z -> max. energy = Z * 100 EeV
# composition = p, He, C, O, Si, Fe with relative abundances from Allard 2006, DOI: 10.1088/1475-7516/2006/09/005

nS = 5 # number of spectral indices between 2 and 3
nP = 5000 # number of particles per spectral index

d = {1:zeros(nS), 2:zeros(nS), 6:zeros(nS), 8:zeros(nS), 14:zeros(nS), 26:zeros(nS)}

print 'simulating for spectral index'
for i in range(nS):
	beta = 2 + i/float(nS - 1)
	print beta

	composition = SourceComposition(10, 100, -beta)
	composition.add(1, 1, 92000)
	composition.add(4, 2, 13000)
	composition.add(12, 6, 447.4)
	composition.add(16, 8, 526.3)
	composition.add(28, 14, 100)
	composition.add(56, 26, 97)

	state = ParticleState()

	for j in range(nP):
		composition.prepare(state)
		z = state.getChargeNumber()
		d[z][i] += 1

	norm = float(d[1][i])
	for z in d.keys():
		d[z][i] /= norm

figure()
beta = linspace(2, 3, nS)
elements = {1:'H', 2:'He', 6:'C', 8:'O', 14:'Si', 26:'Fe'}
for z in d.keys():
	plot(beta, d[z], label=elements[z])

legend(loc = 'lower right')
xlabel(r'Source Spectral Index $\beta$')
ylabel('Relative Integrated Abundance')
xlim(2, 3)
semilogy()
grid()
savefig('SourceCompostion.png')
show()
