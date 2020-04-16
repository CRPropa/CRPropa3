# CRPropa example
# Plots a slice of the JF12 Galactic magnetic field model

from crpropa import *
from pylab import *
from matplotlib.colors import LogNorm


# Components of JF12 GMF model
bField = JF12Field()
bField.randomStriated()
bField.randomTurbulent()

z = 0 # z position [kpc] of slice to plot
N = 241 # samples per direction
samples = (linspace(-20, 20, N, endpoint=True))

B = zeros((N,N))
X, Y = meshgrid(samples, samples)

for i, y in enumerate(samples):
  for j, x in enumerate(samples):
    pos = Vector3d(x, y, z) * kpc
    B[i, j] = bField.getField(pos).getR() / gauss # B in [G]

figure()
maB = ma.masked_array(B, B==0)
pc = pcolor(X, Y, maB, norm=LogNorm(), vmin=1e-9, vmax=1e-4)
gca().set_aspect(1)
cbar = colorbar(pc, shrink=0.8)
cbar.set_label(r'$\vert \vec{B} \vert$ [G]')
plot(-8.5, 0, 'wo')
xlabel('x [kpc]')
ylabel('y [kpc]')
xlim(-20, 20)
ylim(-20, 20)
savefig('JF12.png', bbox_inches='tight')

show()
