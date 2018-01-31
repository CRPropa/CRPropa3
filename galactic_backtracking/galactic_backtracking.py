
# coding: utf-8

# # Galactic backtracking
# The following setup shows how to use CRPropa for backtracking simulations.  
# In the JF12 model the Galaxy is a sphere of 20 kpc radius.
# For the magnetic field we are going to consider the regular component of the JF2012 model. The large-scale (striated) and small-scale (turbulent) random components can optionally be activated with the outcommented sections and a random seed can be set for reproducability.

# In[1]:

from crpropa import *

# magnetic field setup
B = JF12Field()
#seed = 691342
#B.randomStriated(seed)
#B.randomTurbulent(seed)

# simulation setup
sim = ModuleList()
sim.add(PropagationCK(B, 1e-4, 0.1 * parsec, 100 * parsec))
obs = Observer()
obs.add(ObserverLargeSphere(Vector3d(0), 20 * kpc))
# obs.onDetection(TextOutput('galactic_backtracking.txt', Output.Event3D))
sim.add(obs)
print sim


# ## Backtracking a single cosmic ray
# 
# Let's assume we observed a 10 EeV cosmic ray coming from the direction given by longitude and colatitude (1.95, 0.96) radian and want to investigate its direction before having traversed the Galaxy.
# 
# Backtracking corresponds to forward-tracking a particle of the opposite charge, thus we select an anti-proton, which in the HEP ID numbering scheme is denoted by a negative sign.
# Assuming the cosmic ray was a proton the backtracking turns out as follows.

# In[2]:

pid = - nucleusId(1,1)  # (anti-)proton
energy = 10 * EeV
position = Vector3d(-8.5, 0, 0) * kpc
lat = 0.96
lon = 1.95
direction = Vector3d()
direction.setRThetaPhi(1, lat, lon)
p = ParticleState(pid, energy, position, direction)
c = Candidate(p)

sim.run(c)
print c

d1 = c.current.getDirection()  # direction at galactic border
print 'galactic deflection %.2f radian' % direction.getAngleTo(d1)


# ## Backtracking including uncertainties
# The impact of the cosmic ray uncertainties backtracked directions can be investigated with a MC approach. In the following, the cosmic ray energy and direction are varied within the statistical uncertainties before backtracking.

# In[3]:

R = Random()  # CRPropa random number generator

pid = - nucleusId(1,1)
meanEnergy = 10 * EeV
sigmaEnergy = 0.1 * meanEnergy  # 10% energy uncertainty
position = Vector3d(-8.5, 0, 0) * kpc
lat0 = 0.96
lon0 = 1.95
meanDir = Vector3d()
meanDir.setRThetaPhi(1, lat0, lon0)
sigmaDir = 0.002  # 1 degree directional uncertainty

lons, lats = [], []
for i in range(100):
    energy = R.randNorm(meanEnergy, sigmaEnergy)
    direction = R.randVectorAroundMean(meanDir, sigmaDir)

    c = Candidate(ParticleState(pid, energy, position, direction))
    sim.run(c)

    d1 = c.current.getDirection()
    lons.append(d1.getPhi())
    lats.append(d1.getTheta())


# ## (Optional) Plotting
# Finally we are plotting a skymap of the observed direction along with the distribution of directions at the galactic border.

# In[4]:

get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
import numpy as np

# Angle definitions:
# CRPropa uses
#   longitude (phi) [-pi, pi] with 0 pointing in x-direction
#   colatitude (theta) [0, pi] with 0 pointing in z-direction
# matplotlib expects
#   longitude [-pi, pi] with 0 = 0 degrees
#   latitude [pi/2, -pi/2] with pi/2 = 90 degrees (north)
lat0 = np.pi/2 - lat0
lats = np.pi/2 - np.array(lats)

plt.figure(figsize=(12,7))
plt.subplot(111, projection = 'hammer')
plt.scatter(lon0, lat0, marker='+', c='black', s=100)
plt.scatter(lons, lats, marker='o', c='blue', linewidths=0, alpha=0.2)
plt.grid(True)

