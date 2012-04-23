from mpc import *
from pylab import *


def compare(dataFileName, photonField, plotFileName):
    A = genfromtxt(dataFileName, comments="#", unpack=True)
    figure()
    plot(A[0], A[1], label='data')
    
    epp = ElectronPairProduction(photonField)
    
    c = Candidate()
    c.current.setId(getNucleusId(1,1))
    c.setCurrentStep(1 * Mpc)
    
    N = 50
    E = logspace(16,23,N)
    dE = zeros(N)
    
    for i in range(N):
        c.current.setEnergy(E[i] * eV)
        epp.process(c)
        dE[i] = (E[i] - c.current.getEnergy() / eV)
    
    plot(E, dE,'k+', label='simulated', linewidth=2, markeredgewidth=2)
    
    xlabel('Energy [eV]')
    ylabel('Energy Loss Rate [eV / Mpc]')
    legend(loc='lower right')
    grid()
    loglog()
    
    savefig(plotFileName, bbox_inches='tight')


compare(getDataPath("/ElectronPairProduction/cmb.txt"), CMB, 'ElectronPairProduction_CMB.png')
compare(getDataPath("/ElectronPairProduction/ir.txt"), IRB, 'ElectronPairProduction_IRB.png')
compare(getDataPath("/ElectronPairProduction/cmbir.txt"), CMB_IRB, 'ElectronPairProduction_CMB_IRB.png')