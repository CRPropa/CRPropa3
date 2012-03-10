from mpc import *
from pylab import *
import ROOT

ppp = PhotoPionProduction(PhotoPionProduction.CMBIR)
c = Candidate()
c.current.setId(getNucleusId(1,1))
c.current.setEnergy(100 * EeV)
c.setCurrentStep(0.0 * Mpc)

cn = ROOT.TCanvas()
h = ROOT.TH1F('','',50,0,100)

for i in range(10000):
    c.setNextStep(1000 * Mpc)
    c.clearInteractionStates()
    
    ppp.process(c)
    h.Fill(c.getNextStep() / Mpc)

f = ROOT.TF1('f1','expo')
h.Fit(f)
#ROOT.gStyle.setOptFit(1)
#h.Draw()
h.Draw()
cn.SaveAs("PhotoPionProduction_p_100EeV.pdf")

#xlabel('Distance [Mpc]')
#ylabel('Events')
#legend()

E,P,N = genfromtxt(getDataPath('/PhotoPionProduction/cmbir.txt'), unpack=True)
i = argmin(abs(E - 100))

s = 1/-f.GetParameter(1)
ds = f.GetParError(1) * s**2
print s, ds
print 1 / ( P[i-1] + (P[i]-P[i-1]) * (100-E[i-1]) / (E[i]-E[i-1]) )

#savefig('PhotoPionProduction_p_100EeV.png',bbox_inches='tight')

proton, neutron * energien