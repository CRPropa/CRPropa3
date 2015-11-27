# CRPropa test script
# Plots a sclice of the Pshirkov field
#
from crpropa import *
from pylab import *
from matplotlib.colors import LogNorm


# Components of Pshirkov model
# To plot the ASS disk variant of Pshirkov uncomment the corresponding line
bField = PshirkovField()
bField.setUseBSS(True)
#bField.setUseASS(True)

# Check which variant is used, to change the filename
if bField.isUsingASS():
  field = "ASS"
elif bField.isUsingBSS():
  field = "BSS"


z = 0 # z position [kpc] of slice to plot
N = 241 # samples per direction
samples = (linspace(-20, 20, N, endpoint=True))

B = zeros((N,N))
X, Y = meshgrid(samples, samples)

for i,x in enumerate(samples):
  for j,y in enumerate(samples):
    pos = Vector3d(x, y, z) * kpc
    B[i, j] = bField.getField(pos).getR() / gauss # B in [G]

figure()
maB = ma.masked_array(B, B==0)
pc = pcolor(X, Y, maB, norm=LogNorm(), vmin=1e-9, vmax=1e-4)
gca().set_aspect(1)
cbar = colorbar(pc, shrink=0.8)
cbar.set_label(r'$\vert \vec{B} \vert$ [G]')
plot(-8.5, 0, 'wo')
xlabel('x [kpx]')
ylabel('y [kpx]')
xlim(-20, 20)
ylim(-20, 20)
savefig('PT11_%s.png' % field, bbox_inches='tight')

show()
