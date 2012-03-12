import mpc
import pylab as lab
import sys
# make sure root is startet in Batch mode
sys.argv.append( '-b-' )
import ROOT

def getSlope(interaction, energy, charge):
    c = mpc.Candidate()
    c.current.setId(mpc.getNucleusId(1, charge))
    c.current.setEnergy(energy)
    c.setCurrentStep(0.0 * mpc.Mpc)

    h = ROOT.TH1F('', '', 50, 0, 100)

    for i in range(100000):
        c.setNextStep(1000 * mpc.Mpc)
        c.clearInteractionStates()
    
        ppp.process(c)
        h.Fill(c.getNextStep() / mpc.Mpc)

    f = ROOT.TF1('f1', 'expo')
    h.Fit(f)
    
    s = -f.GetParameter(1)
    ds = f.GetParError(1) * s ** 2
    return [s, ds]
    

ppp = mpc.PhotoPionProduction(mpc.PhotoPionProduction.CMBIR)

E, P, N = lab.genfromtxt(mpc.getDataPath('/PhotoPionProduction/cmbir.txt'), unpack=True)
p = lab.zeros(len(E))
n = lab.zeros(len(E))
for i in range(len(E)):
    sds = getSlope(ppp, E[i] * mpc.EeV, 1)
    p[i] = sds[0]
    sds = getSlope(ppp, E[i] * mpc.EeV, 0)
    n[i] = sds[0]

lab.plot(E, n, "b+", label="mpc neutron", markersize=8 )
lab.plot(E, p, "r+", label="mpc proton", markersize=8 )
lab.plot(E, N, "b", label="input neutron"  )
lab.plot(E, P, "r", label="input proton"  )
lab.xlabel('energy [EeV]')
lab.ylabel('rate [1/Mpc]')
lab.legend(loc='lower right')
lab.grid()
lab.semilogx()
lab.savefig('PhotoPionProduction.png',bbox_inches='tight')
#lab.show()
#proton, neutron * energien
