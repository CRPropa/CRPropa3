from mpc import *
from pylab import *


# parameters
lo = 80
hi = 160
origin = Vector3d(lo * Mpc)
N = 128
size = (hi - lo) * Mpc
spacing = size / N
lMin = 2 * spacing
lMax = 16 * spacing
Brms = 1
alpha = -11./3

# create field
field = SPHTurbulentMagneticField(origin, size, N)
print 'initialize'
field.initialize(lMin, lMax, Brms, alpha)
print 'modulate'
field.modulate(getDataPath("SPH/mhd_z.db"), 2./3)
print 'evaluate'

# plot a slice
A = zeros((N,N))
for ix in range(N):
    for iy in range(N):
        b = field.getField(origin + Vector3d(ix, iy, N/2) * spacing)
        A[ix,iy] = b.getMag()

figure()
im = imshow(log10(A), origin='lower', extent=[lo, hi, lo, hi], vmin=-2, vmax=2)
cbar = colorbar(im)
cbar.set_label('$\log_{10}(B/B_{rms})$')
xlabel('x [Mpc]')
ylabel('y [Mpc]')
savefig('SPHTurbulentField_slice.png', bbox_inches='tight')

# sample the effective rms field strengths
brms = 0
r = Random()
for i in range(1000):
    x = Vector3d(r.rand(), r.rand(), r.rand()) * size + origin
    brms += field.getField(x).getMag2()
brms = (brms / 1000)**.5
print 'effective Brms =', brms, 'of nominal Brms'

# sample the coherence length
def getCorrelationCurve():
	x = Vector3d(r.rand(), r.rand(), r.rand()) * size + origin
	d = r.randUnitVectorOnSphere()
	b = field.getField(x)
	b.normalize()
	c = zeros(nSteps)
	for i in range(nSteps):
		b2 = field.getField(x + Vector3d(d.x, d.y, d.z) * i * step)
		b2.normalize()
		c[i] = b.dot(b2)
	return c

nSteps = 100
step = lMax / 2 / nSteps # choose so that lc_unmodulated < nStep * step < size

cMean = zeros(nSteps)
for i in range(1000):
	cMean += getCorrelationCurve() / 1000
lc = (2 * sum(cMean[1:]) + cMean[0]) * step / Mpc

plot(arange(0, nSteps) * step / Mpc, cMean)
xlabel('Linear Separation [Mpc]')
ylabel('Autocorrelation')
text(0.35, 0.95, 'Coherence Length %.1f Mpc'%lc, ha='left', va='top', transform=gca().transAxes)
savefig('SPHTurbulentField_coherenceLength.png', bbox_inches='tight')

