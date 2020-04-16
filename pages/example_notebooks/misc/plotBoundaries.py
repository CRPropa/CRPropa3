# CRPropa example
# Tests and visualizes how periodic and reflective boundaries work.

from crpropa import *
from pylab import *
from matplotlib.patches import Rectangle

p = ParticleState()
p.setPosition(Vector3d(-5, -5, 0))
p.setDirection(Vector3d(1,0.25,0))

N = 45
x, y, sx, sy = zeros((4,N))

# Periodic boundaries
m = ModuleList()
m.add(SimplePropagation(0, 1))
m.add(PeriodicBox(Vector3d(-10), Vector3d(20)))

c = Candidate(p)
for i in range(N):
	m.process(c)
	pos = c.current.getPosition()
	x[i], y[i] = pos.x, pos.y
	src = c.source.getPosition()
	sx[i], sy[i] = src.x, src.y

figure(figsize=(8,6))
ax1 = subplot(211, aspect='equal')
ax1.add_patch(Rectangle((-10,-10), 20, 20, facecolor="lightskyblue"))
ax1.plot(x, y, 'b.', ms=5)
ax1.plot(sx, sy, 'r*', markeredgewidth=0, ms=10)
ax1.set_xlim((-51, 31))
ax1.set_ylim((-11, 11))
ax1.grid()
ax1.set_title('Periodic Boundaries')

# Reflective boundaries
m = ModuleList()
m.add(SimplePropagation(0, 1))
m.add(ReflectiveBox(Vector3d(-10), Vector3d(20)))

c = Candidate(p)
for i in range(N):
	m.process(c)
	pos = c.current.getPosition()
	x[i], y[i] = pos.x, pos.y
	src = c.source.getPosition()
	sx[i], sy[i] = src.x, src.y

ax2 = subplot(212, aspect='equal')
ax2.add_patch(Rectangle((-10,-10), 20, 20, facecolor="lightskyblue"))
ax2.plot(x, y, 'b.', ms=5)
ax2.plot(sx, sy, 'r*', markeredgewidth=0, ms=10)
ax2.set_xlim((-51, 31))
ax2.set_ylim((-11, 11))
ax2.grid()
ax2.set_title('Reflective Boundaries')

savefig('Boundaries.png', bbox_inches='tight')
show()
