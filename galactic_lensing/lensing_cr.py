
# coding: utf-8

# # Galactic Lensing of Simulated Cosmic Rays
# 
# Deflection in galactic magnetic fields can be accounted for using galactic lenses. To use the lenses efficiently, the ouput of a CRPropa simulation needs first to be filled in probability maps. These maps are then transformed by the lenses according to the rigidity of the particles. From the maps then finally new particles can be generated. 
# 
# ## Input Results from Extragalactic Simulation

# In[1]:

import crpropa
import numpy as np
# read data from crpropa output into container. 
# The data is weighted with the source energy ~E**-1
M = crpropa.ParticleMapsContainer()
crdata = np.genfromtxt("crpropa_output.txt")
id = crdata[:,0]
E = crdata[:,3] * crpropa.EeV
E0 = crdata[:,4] * crpropa.EeV
px = crdata[:,11]
py = crdata[:,12]
pz = crdata[:,13]
galactic_longitude = np.arctan2(-1. * py, -1. *px)
galactic_latitude = np.pi / 2 - np.arccos( -pz / (px*px + py*py+ pz*pz) )

M.addParticles(id, E, galactic_longitude, galactic_latitude, E0**-1)


# Alternatively, data can be added manually. 
# This provides freedom to adapt to customized weights used in the simulation
for i in xrange(1000):
    particleId = crpropa.nucleusId(1,1)
    energy = 10 * crpropa.EeV
    galCenter = crpropa.Vector3d(-1,0,0)
    momentumVector = crpropa.Random.instance().randFisherVector(galCenter, 200)
    M.addParticle(particleId, energy, momentumVector)   



# ## Plot Maps
# 
# The probability maps for the individual particles can be accessed directly and plotted with healpy. 
# 

# In[2]:

get_ipython().magic('matplotlib inline')

import healpy
from pylab import *


#stack all maps
crMap = zeros(49152)
for pid in M.getParticleIds():
    energies = M.getEnergies(int(pid))
    for i, energy in enumerate(energies):
        crMap += M.getMap(int(pid), energy * crpropa.eV )

#plot maps using healpy
healpy.mollview(map=crMap, title='Unlensed')
savefig('unlensed_map.png')


# ## Apply Galactic Lenses
# 
# To apply a lens to a map, a lens object needs to be created and then applied to the map. The normalization of the lens ensures that the lens does not distort the spectra.

# In[3]:

get_ipython().magic('matplotlib inline')

lens = crpropa.MagneticLens('pathto/lens.cfg')
lens.normalizeLens()
M.applyLens(lens)

#stack all maps
crMap = zeros(49152)
for pid in M.getParticleIds():
    energies = M.getEnergies(int(pid))
    for i, energy in enumerate(energies):
        crMap += M.getMap(int(pid), energy * crpropa.eV )

#plot maps using healpy
healpy.mollview(map=crMap, title='Lensed')
savefig('lensed_map.png')


# ## Generate Particles from Map
# 
# If neeeded, sets of individual particles can then be generated from the maps under preservation of the relative (transformed) input weights. Note that the numbr of particles you can draw from a map depends on the quality of sampling on the input pdf. As a rule of thumb it is most likely not safe to draw more particles from the map than have been used to create the map. 

# In[4]:

pids, energies, lons, lats = M.getRandomParticles(10)

# create a scatter plot of the particles
subplot(111, projection='hammer')
scatter(lons, lats, c=log10(energies), lw=0)
savefig('scatterd_particles.png')

