from mpc import *
from pylab import *

D = zeros((27, 31))
D2 = zeros((27, 31))

for z in range(0, 27):
    for n in range(0, 31):
        try:
            D[z, n] = getNucleusMass(getNucleusId(n + z, z)) / amu
            D2[z, n] = D[z, n] - (z + n) - z * mass_electron / amu
        except:
            continue

# nuclear masses
figure()
im = imshow(ma.masked_array(D, D==0), aspect='equal', interpolation='nearest', origin='lower')
cbar = colorbar(im)
cbar.set_label('Mass [amu]')
xlabel('Neutrons')
ylabel('Protons')
grid()
savefig('NuclearMass_mass.png', bbox_inches='tight')
