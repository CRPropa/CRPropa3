
# coding: utf-8

# # Photon Propagation
# 
# There are two main ways to propagate electromagenetic particle (EM particles: photons, electrons, positrons) in CRPropa.
# 1) propagation as part of the CRPropa simulation chain
# 2) propagation outside of the CRPropa simulation chain
# 
# The following describes option 2, for which CRPropa provides three functions.
# EM particles can be either propagated individually using the external EleCa code (suitable for high energies), or their spectra can be propapagated with the transport code DINT (suitable for low energies).
# Alternatively a combined option is available that processes high energy photons with Eleca and the calculates the resulting spectra with DINT down to low energies.
# 
# All three functions take as input a plain-text file with EM particles in `PhotonOutput1D` format.
# In the following examples the input file "photon_monoenergetic_source.dat" contains 1000 photons with E = 50 EeV from a photon source in 4 Mpc distance.
# 
# The last example "Photons from Proton Propagation" shows how to obtain secondary EM particles from a simulation of hadronic cosmic rays.
# 
# Note that the differing results in EleCa (and correspondingly the high energy part of the combined option) are due to an incorrect sampling of the background photon energies in EleCa. The EleCa support will be removed in the near future.
# 
# 

# ### Propagation with Eleca
# 

# In[1]:

import crpropa

# Signature: ElecaPropagation(inputfile, outputfile, showProgress=True, lowerEnergyThreshold=5*EeV, magneticFieldStrength=1*nG, background="ALL")
crpropa.ElecaPropagation("photon_monoenergetic_source.dat", "photons_eleca.dat", True, 0.1*crpropa.EeV, 0.1*crpropa.nG)



# ### Propagation with DINT
# 

# In[2]:

import crpropa

# Signature: DintPropagation(inputfile, outputfile, IRFlag=4, RadioFlag=4, magneticFieldStrength=1*nG, aCutcascade_Magfield=0)
crpropa.DintPropagation("photon_monoenergetic_source.dat", "spectrum_dint.dat", 4, 4, 0.1*crpropa.nG)



# ### Combined Propagation

# In[3]:

import crpropa

# Signature: DintElecaPropagation(inputfile, outputfile, showProgress=True, crossOverEnergy=0.5*EeV, magneticFieldStrength=1*nG, aCutcascade_Magfield=0)
crpropa.DintElecaPropagation("photon_monoenergetic_source.dat", "spectrum_dint_eleca.dat", True, 0.5*crpropa.EeV, 0.1*crpropa.nG)


# ### (Optional) Plotting of Results

# In[6]:

get_ipython().magic('matplotlib inline')
from pylab import *

figure(figsize=(6,6))

loglog(clip_on=False)
yscale("log", nonposy='clip')
xlabel('Energy [eV]')
ylabel ('$dN/dE\; E^2$  [a.u.]')

# Plot the Eleca spectrum
elecaPhotons = genfromtxt("photons_eleca.dat")
binEdges = 10**arange(12, 24, .1)
logBinCenters = log10(binEdges[:-1]) + 0.5 * (log10(binEdges[1:]) - log10(binEdges[:-1]))
binWidths = (binEdges[1:] - binEdges[:-1])
data = histogram(elecaPhotons[:,1] * 1E18, bins=binEdges)
J = data[0] / binWidths
dJ = J / sqrt(data[0])
E = 10**logBinCenters
step(E, J * E**2,  c='m', label='EleCa')

#Plot the dint spectrum
data = genfromtxt("spectrum_dint.dat", names=True)
lE = data['logE']
E  = 10**lE
dE = 10**(lE + 0.05) - 10**(lE - 0.05)
J  = data['photons'] / dE
step(E, J * E**2 , c='b', ls='-', where='mid', label='DINT')

#Plot the combined dint+Eleca spectrum
data = genfromtxt("spectrum_dint_eleca.dat", names=True)
lE = data['logE']
E  = 10**lE
dE = 10**(lE + 0.05) - 10**(lE - 0.05)
J  = data['photons'] / dE
step(E, J * E**2 , c='r', ls='-', where='mid', label='Combined')

# Nice limits
xlim(1e14, 1E21)
ylim(ymin=1E17)
legend(loc='upper left')





# ## Photons from Proton Propagation
# 
# The generaton of photons has to be enabled for the individual energy-loss processes in the module chain. Also, a photon output module has to be added:

# In[ ]:

from crpropa import *   

# source setup  
source = Source()       
source.add(SourceParticleType(nucleusId(1, 1))) 
source.add(SourcePowerLawSpectrum(10 * EeV, 100 * EeV, -2))     
source.add(SourceUniform1D(3 * Mpc, 100.00001 * Mpc))   

# setup module list for proton propagation
m = ModuleList()
m.add(SimplePropagation(0, 10 * Mpc))
m.add(MinimumEnergy(1 * EeV))

# observer
obs = Observer()
obs.add( ObserverPoint() )
obs.onDetection( TextOutput('proton_output.txt', Output.Event1D) )
m.add(obs)

# secondary electrons are disabled here
m.add(ElectronPairProduction(CMB, False))        
# enable secondary photons
m.add(PhotoPionProduction(CMB, True))   

# add a photon output
m.add(PhotonOutput1D("photon_output.txt"))

# run simulation        
m.run(source, 10000, True)


# The file 'photon_output' will contain approximately 600 secondary particles, half of them photons and can be processed as the photon example above.
#ID     E       D       pID     pE      iID     iE
#
# ID          Id of particle (photon, electron, positron)
# E           Energy [EeV]
# D           Comoving trajectory length [Mpc]
# pID         Id of parent particle
# pE          Energy [EeV] of parent particle
# iID         Id of source particle
# iE          Energy [EeV] of source particle
#
 -11    6.77482  67.1766        1000000010       75.9465        1000010010       91.8829
  11    0.0191816        20.3029        1000010010       65.2432        1000010010       91.8829
 -11    0.346937         30.1839        1000000010       61.6767        1000010010       85.4494
 -11    0.414738         20.0958        1000000010       54.7093        1000010010       63.1855
 -11    0.967904          9.0220        1000000010       81.2254        1000010010       94.8547
  22    14.0028  50.9686        1000010010       72.0608        1000010010       88.8126
  22    8.24537  43.0673        1000010010       70.7568        1000010010       92.8280
  22    0.771388         50.9686        1000010010       72.0608        1000010010       88.8126
  22    10.5142  43.0673        1000010010       70.7568        1000010010       92.8280
 -11    2.61693  50.9686        1000000010       64.1018        1000010010       88.8126
 -11    2.0019   33.0673        1000000010       60.4899        1000010010       92.8280
  22    14.31     3.0673        1000000010       39.5666        1000010010       92.8280
  22    6.61331   3.0673        1000000010       39.5666        1000010010       92.8280
 -11    2.18491  52.5561        1000000010       53.0487        1000010010       63.5297
 -11    4.29165  33.0320        1000000010       73.2522        1000010010       83.5731
 -11    0.441877         62.9243        1000000010       90.4563        1000010010       99.4810
  11    1.45     46.9630        1000010010       83.2926        1000010010       99.4810
  22    0.876782         36.3717        1000010010       55.3764        1000010010       66.9174
  22    10.1596  36.3717        1000010010       55.3764        1000010010       66.9174
 -11    0.456958         64.0524        1000000010       55.8222        1000010010       64.8143
 -11    2.64216  37.4412        1000000010       74.5617        1000010010       90.1071
 -11    5.70889  38.7002        1000000010       49.5142        1000010010       64.7278
 -11    1.05272  38.6030        1000000010       47.4831        1000010010       59.5987
...
# In[ ]:



