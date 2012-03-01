import sys
sys.path.append("/home/walz/software/mpc/build")

from mpc import *

c = Candidate()
c.current.setId(1001001000)
c.current.setEnergy(100 * EeV)
c.setCurrentStep(1 * Mpc)

epp = ElectronPairProduction(ElectronPairProduction.CMBIR)
epp.process(c)
