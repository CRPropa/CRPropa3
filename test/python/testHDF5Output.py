# CRPropa test script
# Tests and visualizes how periodic and reflective boundaries work.
#
from crpropa import *


p = ParticleState()
p.setPosition(Vector3d(-5, -5, 0) * Mpc)
p.setDirection(Vector3d(1,0.25,0))

N = 4500000

# Periodic boundaries
m = ModuleList()
m.add(SimplePropagation(0, 1))
m.add(PeriodicBox(Vector3d(-10), Vector3d(20)))
m.add(HDF5Output("test.h5"))

print ("Run...") 
c = Candidate(p)
for i in range(N):
	if i % 1000 == 0:
		print (i)
	m.process(c)
