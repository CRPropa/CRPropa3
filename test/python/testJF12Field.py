from mpc import *
from pylab import *
from matplotlib.colors import LogNorm


bField = JF2012Field()
#bField.randomStriated()
#bField.randomTurbulent()

N = 241
lx = (linspace(-20, 20, N, endpoint=True))
z = 0
B, X, Y = zeros((3, N,N))

for ix in range(N):
  for iy in range(N):
    b = bField.getField(Vector3d(lx[ix], lx[iy], z) * kpc)
    B[ix, iy] = b.getMag() / gauss # B in [G]
    X[ix, iy] = lx[ix]
    Y[ix, iy] = lx[iy]

figure()
maB = ma.masked_array(B, B==0)
pc = pcolor(X, Y, maB, norm=LogNorm(), vmin=1e-9, vmax=1e-4)
gca().set_aspect(1)
cbar = colorbar(pc, shrink=0.8)
cbar.set_label(r'$\vert \vec{B} \vert$ [G]')
plot(-8.5, 0, 'wo')
xlabel('x [kpx]')
ylabel('y [kpx]')
show()
savefig('JF12.png', bbox_inches='tight')
