from mpc import *
from pylab import *

def getRate(module, energy, charge):
    c = Candidate()
    c.current.setId(nucleusId(1, charge))
    c.current.setEnergy(energy)
    
    N = 5000
    l = 0
    for i in range(N):
        s = InteractionState()
        module.setNextInteraction(c, s)
        l += s.distance / Mpc / N
    
    return 1 / l
    
def compare(dataFileName, photonField, plotFileName):
    E_data, P_data, N_data = genfromtxt(dataFileName, unpack=1)
    E = E_data[1:-1:5]
    p = zeros(len(E))
    n = zeros(len(E))

    ppp = PhotoPionProduction(photonField)
    for i,energy in enumerate(E*EeV):
        p[i] = getRate(ppp, energy, 1)
        n[i] = getRate(ppp, energy, 0)

    figure()
    plot(E_data, P_data, "r", label="Proton data table")
    plot(E, p, 'k+', label="Proton")
    plot(E_data, N_data, "b", label="Neutron data table")
    plot(E, n, 'k.', label="Neutron")
    xlabel('Energy [EeV]')
    ylabel('Rate [1/Mpc]')
    legend(loc='center right')
    grid()
    loglog()
    xlim(10, 1e5)
    ylim(1e-4, 1)
    savefig(plotFileName)

compare(getDataPath('photopion_CMB.txt'), CMB, 'PhotoPionProduction_CMB.png')
compare(getDataPath('photopion_IRB.txt'), IRB, 'PhotoPionProduction_IRB.png')

show()
