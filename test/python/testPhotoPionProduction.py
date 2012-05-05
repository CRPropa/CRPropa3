from mpc import *
from pylab import *

def getRate(module, energy, charge):
    c = Candidate()
    c.current.setId(getNucleusId(1, charge))
    c.current.setEnergy(energy)
    
    N = 5000
    l = 0
    for i in range(N):
        s = InteractionState()
        module.setNextInteraction(c, s)
        l += s.distance / Mpc / N
    
    return 1 / l
    
def compare(type, name):
    print "compare ", name
    ppp = PhotoPionProduction(type)
    E_data, P_data, N_data = genfromtxt(getDataPath('/PhotoPionProduction/PPtable_' + name + '.txt'), unpack=True)

    E = E_data[1:-1:5]
    p = zeros(len(E))
    n = zeros(len(E))
    for i,energy in enumerate(E*EeV):
        p[i] = getRate(ppp, energy, 1)
        n[i] = getRate(ppp, energy, 0)
    figure()
    plot(E_data, P_data, "r", label="Proton Data")
    plot(E, p, 'k+', label="Proton Simulated")
    plot(E_data, N_data, "b", label="Neutron Data")
    plot(E, n, 'k.', label="Neutron Simulated")
    xlabel('Energy [EeV]')
    ylabel('Rate [1/Mpc]')
    legend(loc='center right')
    grid()
    loglog()
    xlim(10, 1e5)
    ylim(1e-4, 1)
    savefig('PhotoPionProduction_' + name + '.png', bbox_inches='tight')

compare(CMB, "CMB")
compare(IRB, "IRB")
