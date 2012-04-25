from mpc import *
from pylab import *
import ROOT
ROOT.gROOT.SetBatch(True)


def getSlope(module, energy, charge, name):
    c = Candidate()
    c.current.setId(getNucleusId(1, charge))
    c.current.setEnergy(energy)
    
    l = []
    for i in range(5000):
        s = InteractionState()
        module.setNextInteraction(c, s)
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

    E = E_data[1:-1:5]
    p = zeros(len(E))
    n = zeros(len(E))
    for i,energy in enumerate(E*EeV):
        p[i] = getSlope(ppp, energy, 1, name)
        n[i] = getSlope(ppp, energy, 0, name)
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
