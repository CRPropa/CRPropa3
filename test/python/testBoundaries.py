import pylab as pl
import matplotlib.patches as mpatches
from mpc import *


p = ParticleState()
p.setPosition(Vector3d(-5, -5, 0))
p.setDirection(Vector3d(1,0.25,0))

N = 45
x, y, sx, sy = pl.zeros((4,N))


# simulate periodic boundaries
m1 = ModuleList()
m1.add(SimplePropagation(0, 1))
m1.add(PeriodicBox(Vector3d(-10), Vector3d(20)))

c = Candidate(p)
for i in range(N):
	m1.process(c)
	pos = c.current.getPosition()
	x[i], y[i] = pos.x, pos.y
	src = c.initial.getPosition()
	sx[i], sy[i] = src.x, src.y

pl.figure(figsize=(8,6))
ax1 = pl.subplot(211, aspect='equal')
ax1.add_patch(mpatches.Rectangle((-10,-10), 20, 20, facecolor="lightskyblue"))
ax1.plot(x, y, 'b.', ms=5)
ax1.plot(sx, sy, 'r*', markeredgewidth=0, ms=10)
ax1.set_xlim((-51, 31))
ax1.set_ylim((-11, 11))
ax1.grid()
ax1.set_title('Periodic Boundaries')


# simulate reflective boundaries
m2 = ModuleList()
m2.add(SimplePropagation(0, 1))
m2.add(ReflectiveBox(Vector3d(-10), Vector3d(20)))

c = Candidate(p)

for i in range(N):
	m2.process(c)
	pos = c.current.getPosition()
	x[i], y[i] = pos.x, pos.y
	src = c.initial.getPosition()
	sx[i], sy[i] = src.x, src.y

ax2 = pl.subplot(212, aspect='equal')
ax2.add_patch(mpatches.Rectangle((-10,-10), 20, 20, facecolor="lightskyblue"))
ax2.plot(x, y, 'b.', ms=5)
ax2.plot(sx, sy, 'r*', markeredgewidth=0, ms=10)
ax2.set_xlim((-51, 31))
ax2.set_ylim((-11, 11))
ax2.grid()
ax2.set_title('Reflective Boundaries')

pl.savefig('Boundaries.png', bbox_inches='tight')
pl.show()
