# CRPRopa test script
# Shows how to use the energy loss length method, provided by the interaction modules
#
from crpropa import *
from pylab import *

pd1 = PhotoDisintegration(CMB)
pd2 = PhotoDisintegration(IRB)
pp1 = PhotoPionProduction(CMB)
pp2 = PhotoPionProduction(IRB)
ep1 = ElectronPairProduction(CMB)
ep2 = ElectronPairProduction(IRB)

def parse_pid(pid):
    return 'Z=%i, A=%i'%((pid//10000)%1000, (pid%10000)/10)

def invAdd(x1, x2):
    return 1. / (1. / x1 + 1. / x2)

pid = nucleusId(4, 2)
mass = nucleusMass(pid)
energies = logspace(1, 4, 50)
L = zeros((3, 50))
z = 0.

for i, E in enumerate(energies * EeV):
    L[0,i] = invAdd(pd1.lossLength(pid, E, z), pd2.lossLength(pid, E, z))
    L[1,i] = invAdd(pp1.lossLength(pid, E, z), pp2.lossLength(pid, E, z))
    lorentzfactor = E / (mass * c_squared)
    L[2,i] = invAdd(ep1.lossLength(pid, lorentzfactor, z), ep2.lossLength(pid, lorentzfactor, z))
L /= Mpc


figure()
plot(energies, L[0], label='PDis')
plot(energies, L[1], label='PPion')
plot(energies, L[2], label='EPair')
legend(frameon=0, loc='lower left')
text(0.05, 0.25, parse_pid(pid), transform=gca().transAxes)
xlabel('Energy [EeV]')
ylabel('Energy Loss Length [Mpc]')
ylim(1, 1e4)
loglog()
savefig('EnergyLossLength_%i.png'%pid, bbox_inches='tight')
show()
