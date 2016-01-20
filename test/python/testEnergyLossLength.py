# CRPRopa test script
# Shows how to use the energy loss length method, provided by the interaction modules
#
from crpropa import *
from pylab import *

z = 0
pid = nucleusId(4, 2)
mass = nuclearMass(pid)
energies = logspace(1, 4, 50)


pd1 = PhotoDisintegration(CMB)
pd2 = PhotoDisintegration(IRB)
pp1 = PhotoPionProduction(CMB)
pp2 = PhotoPionProduction(IRB)
ep1 = ElectronPairProduction(CMB)
ep2 = ElectronPairProduction(IRB)

pp1Loss, pp2Loss = zeros(50), zeros(50)
pd1Loss, pd2Loss = zeros(50), zeros(50)
ep1Loss, ep2Loss = zeros(50), zeros(50)

for i, E in enumerate(energies * EeV):
    lf = E / (mass * c_squared)  # Lorentz factor

    # pion production
    pp1Loss[i] = pp1.lossLength(pid, lf, z)
    pp2Loss[i] = pp2.lossLength(pid, lf, z)
    # disintegration
    pd1Loss[i] = pd1.lossLength(pid, lf, z)
    pd2Loss[i] = pd2.lossLength(pid, lf, z)
    # pair production
    ep1Loss[i] = ep1.lossLength(pid, lf, z)
    ep2Loss[i] = ep2.lossLength(pid, lf, z)

# inversely add loss lengths
ppLoss = 1/(1/pp1Loss + 1/pp2Loss) / Mpc
pdLoss = 1/(1/pd1Loss + 1/pd2Loss) / Mpc
epLoss = 1/(1/ep1Loss + 1/ep2Loss) / Mpc


figure()
plot(energies, pdLoss, label='PDis')
plot(energies, ppLoss, label='PPion')
plot(energies, epLoss, label='EPair')
legend(frameon=0, loc='lower left')
text(0.05, 0.25,
    'Z=%i, A=%i'%(chargeNumber(pid), massNumber(pid)),
    transform=gca().transAxes)
xlabel('Energy [EeV]')
ylabel('Energy Loss Length [Mpc]')
ylim(1, 1e4)
loglog()
savefig('EnergyLossLength_%i.png'%pid, bbox_inches='tight')
show()
