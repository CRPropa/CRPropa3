import os
dataPath = "/home/walz/software/mpc/data"
os.environ["MPC_DATA_PATH"] = dataPath

import sys
sys.path.append("/home/walz/software/mpc/build")
from mpc import *

from pylab import *
import ROOT



ppp = PhotoPionProduction(PhotoPionProduction.CMBIR)
c = Candidate()
c.current.setId(getNucleusId(1,1))
c.current.setEnergy(100 * EeV)
c.setCurrentStep(0 * Mpc)

h = ROOT.TH1F('','',50,0,400)

for i in range(10000):
    c.setNextStep(1000 * Mpc)
    c.clearInteractionStates()
    
    ppp.process(c)
    h.Fill(c.getNextStep() / Mpc)

f = ROOT.TF1('f1','expo')
h.Fit(f)
s = 1/-f.GetParameter(1)
ds = f.GetParError(1) * s**2

#xlabel('Distance [Mpc]')
#ylabel('Events')
#legend()

E,P,N = genfromtxt(dataPath+'/PhotoPionProduction/cmbir.txt', unpack=True)
i = argmin(abs(E - 100))

print s, ds
print 1 / ( P[i-1] + (P[i]-P[i-1]) * (100-E[i-1]) / (E[i]-E[i-1]) )

#savefig('PhotoPionProduction_p_100EeV.png',bbox_inches='tight')
