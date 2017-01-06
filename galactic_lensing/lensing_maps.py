
# coding: utf-8

# # Galactic Lensing of Maps

# In[1]:

get_ipython().magic('matplotlib inline')

import crpropa
import healpy
from pylab import *
import copy

lens = crpropa.MagneticLens('pathto/lens.cfg')


# In[15]:

# generate random input map
l = [1,1,0,0]
inputMap = abs(healpy.synfast(l, nside=64))

#plot map
healpy.mollview(map=inputMap, title='Unlensed')
savefig('lensed_map.png')


# In[16]:

#copy input data as we need it later
outputMap = copy.copy(inputMap)

#transform map
lens.transformModelVector(outputMap, 5 * crpropa.EeV)

#plot transformed map
healpy.mollview(map=outputMap, title='Lensed')


# In[17]:

#Calculate power spectra
lmax = 20
mmUnlensed = healpy.anafast(inputMap, lmax=lmax)
mmLensed = healpy.anafast(outputMap, lmax=lmax)

#plot power spectra
figure()
l = arange(0, 20.01)
step(l, mmUnlensed / mmUnlensed[0], c='k', label='Unlensed', where='mid')
step(l, mmLensed / mmLensed[0], c='r', label='Lensed', where='mid')
semilogy()
ylim(1E-7,1)
xlabel('$l$')
ylabel(' $C_l / C_0$')
legend()


# In[17]:



