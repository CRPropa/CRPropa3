'''
Simulate the p/n-ratio of the final baryon states in SOPHIA photo-pion interactions
'''
from mpc import *
from pylab import *

spp = SophiaPhotoPionProduction()
c = Candidate()
s = InteractionState()

def getPNRatio(E):
	nP, nN = 0, 0
	for i in range(2000):
		c.current.setId(getNucleusId(1, 1))
		c.current.setEnergy(E * EeV)
		spp.setNextInteraction(c, s)
		spp.performInteraction(c)
		if c.current.getId() == getNucleusId(1,1):
			nP += 1.
		if c.current.getId() == getNucleusId(1,0):
			nN += 1.
	return nP / nN

E = logspace(log10(50), 4, 20)
Rpn = zeros(20)
for i in range(20):
	print i
	Rpn[i] = getPNRatio(E[i])

plot(E, Rpn)
grid()
xlabel('Energy of Incident Proton [EeV]')
ylabel('Proton / Neutron Ratio')
semilogx()
savefig('PhotoPionProduction_pnRatio.png', bbox_inches='tight')
