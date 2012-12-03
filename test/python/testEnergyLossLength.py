from pylab import *
from mpc import *

pd1 = PhotoDisintegration(CMB)
pd2 = PhotoDisintegration(IRB)
pp1 = PhotoPionProduction(CMB)
pp2 = PhotoPionProduction(IRB)
ep = ElectronPairProduction(CMB_IRB)

def pdLoss(pid, E):
  l1 = pd1.energyLossLength(pid, E)
  l2 = pd2.energyLossLength(pid, E)
  return 1/(1/l1 + 1/l2)

def ppLoss(pid, E):
  l1 = pp1.energyLossLength(pid, E)
  l2 = pp2.energyLossLength(pid, E)
  l = 1/(1./l1 + 1./l2)
  if l > 1e4 * Mpc:
    return nan
  return l

def parse_pid(pid):
  return 'Z=%i, A=%i'%((pid//10000)%1000, (pid%10000)/10)



pid = nucleusId(56, 26)

E = logspace(1, 4) * EeV
L = zeros((3, 50))

for i, energy in enumerate(E):
  L[0,i] = pdLoss(pid, energy)
  L[1,i] = ppLoss(pid, energy)
  L[2,i] = ep.energyLossLength(pid, energy)

E /= EeV
L /= Mpc



figure()
plot(E, L[0], label='PDis')
plot(E, L[1], label='PPion')
plot(E, L[2], label='EPair')
legend(frameon=0, loc='lower left')
text(0.05, 0.25, parse_pid(pid), transform=gca().transAxes)
xlabel('Energy [EeV]')
ylabel('Energy Loss Length [Mpc]')
ylim(0.1, 1e4)
loglog()
savefig('EnergyLossLength_%i.png'%pid)
show()
