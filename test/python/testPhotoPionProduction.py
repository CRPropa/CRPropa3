from mpc import *
from pylab import *
import ROOT

def getSlope(interaction, energy, charge, name):
    c = Candidate()
    c.current.setId(getNucleusId(1, charge))
    c.current.setEnergy(energy)
    s = InteractionState()

    l = []
    for i in range(10000):
        c.clearInteractionStates()
        interaction.process(c)
        c.getInteractionState('PhotoPionProduction:'+name, s)
        l.append(s.distance / Mpc)

    h = ROOT.TH1F('', '', 100, 0, max(l))

    for element in l:
        h.Fill(element)

    f = ROOT.TF1('f1', 'expo')
    h.Fit(f, "q")
    s = -f.GetParameter(1)
    ds = f.GetParError(1) * s ** 2
    return s
    
def compare(type, name):
    print "compare ", name
    ppp = PhotoPionProduction(type)
    E_data, P_data, N_data = genfromtxt(getDataPath('/PhotoPionProduction/PPtable_' + name + '.txt'), unpack=True)

    E = logspace(1.5,5,10) * EeV
    p = zeros(len(E))
    n = zeros(len(E))
    for i,energy in enumerate(E*EeV):
        p[i] = getSlope(ppp, energy, 1, name)
        n[i] = getSlope(ppp, energy, 0, name)
    plot(E_data, P_data, "r", label="Proton Data")
    plot(E, p, 'k+', label="Proton Simulated")
    plot(E_data, N_data, "b", label="Neutron Data")
    plot(E, n, 'k.', label="Neutron Simulated")
    xlabel('Energy [EeV]')
    ylabel('Rate [1/Mpc]')
    legend(loc='center right')
    grid()
    loglog()
    ylim(1e-4, 1)
    savefig('PhotoPionProduction_' + name + '.png', bbox_inches='tight')
    close()

compare(CMB, "CMB")
compare(IRB, "IRB")
compare(CMB_IRB, "CMB_IRB")
