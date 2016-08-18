# CRPRopa test script
# Plot and self consistency test for electron pair production.
#
import matplotlib
matplotlib.use('Agg')

from crpropa import *
from pylab import *


def compare(dataFileName, photonField, plotFileName):
    epp = ElectronPairProduction(photonField)
    c = Candidate()
    c.current.setId(nucleusId(1,1))
    c.setCurrentStep(1 * Mpc)

    d = genfromtxt(dataFileName, unpack=True)
    gamma = 10**d[0]  # Lorentz factor
    loss1 = d[1]  # loss length
    loss2 = zeros_like(loss1)

    for i, g in enumerate(gamma):
        c.current.setLorentzFactor(g)
        epp.process(c)
        loss2[i] = (g - c.current.getLorentzFactor()) / g

    figure()
    plot(gamma, loss1, label='Datatable')
    plot(gamma, loss2,'k+', label='Simulated', linewidth=2, markeredgewidth=2)
    xlabel(r'$\log_{10}(\gamma)$')
    ylabel(r'Loss length $\frac{1}{\gamma} d\gamma/dx$ [1/Mpc]')
    legend(loc='lower right')
    grid()
    loglog()
    savefig(plotFileName, bbox_inches='tight')


# CMB
compare(getDataPath('epp_CMB.txt'), CMB, 'ElectronPairProduction_CMB.png')
# Dominguez '11 EBL
compare(getDataPath('epp_IRB_Dominguez11.txt'), IRB_Dominguez11, 'ElectronPairProduction_IRB_Dominguez11.png')
show()
