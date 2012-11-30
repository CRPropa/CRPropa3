from mpc import *
from pylab import *
from matplotlib.colors import LogNorm


bField = JF2012Field()
#bField.randomStriated()
#bField.randomTurbulent()

N = 241
x = (linspace(-20, 20, N, endpoint=True)) * kpc
y = (linspace(-20, 20, N, endpoint=True)) * kpc
z = 0 * kpc

B = zeros((N,N))
for ix in range(N):
  for iy in range(N):
    b = bField.getField(Vector3d(x[ix], y[iy], z))
    B[iy, ix] = b.getMag() / gauss # B in [G], rows = y, columns = x

maB = ma.masked_array(B, B==0)
X, Y = meshgrid(x / kpc, y / kpc)

figure()
pc = pcolor(X, Y, maB, norm=LogNorm(), vmin=1e-9, vmax=1e-4)
gca().set_aspect(1)
cbar = colorbar(pc, shrink=0.8)
cbar.set_label(r'$\vert \vec{B} \vert$ [G]')
xlabel('x [kpx]')
ylabel('y [kpx]')
show()
savefig('JF12.png', bbox_inches='tight')
