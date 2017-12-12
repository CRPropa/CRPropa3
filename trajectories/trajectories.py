
# coding: utf-8

# # 3D trajectories in a turbulent field
# 
# The following simulation tracks a single UHE nucleus and its secondary nucleons/nuclei through a turbulent magnetic field.
# 
# First we create a random realization of a turbulent field with a Kolmogorov power spectrum on 60-800 kpc lengthscales and and an RMS field strength of 8 nG.
# The field is stored on a $256^3$ grid with 30 kpc grid spacing, and thus has an extent of $(256 \cdot 30 \rm{kpc})^3$.
# The field is by default periodically repeated in space to cover an arbitrary volume.
# 
# The chosen grid size consumes only very little memory. For practical purposes a larger grid is advised in order to represent more variations of turbulent modes, provide a larger turbulent range, or a higher resolution.

# In[1]:


from __future__ import print_function
from crpropa import *

randomSeed = 42
vgrid = VectorGrid(Vector3d(0), 256, 30*kpc)
initTurbulence(vgrid, 8*nG, 60*kpc, 800*kpc, -11./3., randomSeed)
Bfield = MagneticFieldGrid(vgrid)

# print some properties of our field
print('Lc = {:.1f} kpc'.format(turbulentCorrelationLength(60, 800, -11./3.))) # correlation length
print('<B^2> = {:.1f} nG'.format((rmsFieldStrength(vgrid) / nG)))   # RMS
print('<|B|> = {:.1f} nG'.format((meanFieldStrength(vgrid) / nG)))  # mean
print('B(10 Mpc, 0, 0) = {}'.format(Bfield.getField(Vector3d(10,0,0) * Mpc) / nG, 'nG'))


# ### Saving and loading fields
# 
# In addition to creating random turbulent fields, we can also load and save custom magnetic field grids.
# As input and output we currently support binary files in single precision and ASCII files.

# In[2]:


# save the field
# format: (Bx, By, Bz)(x, y, z) with z changing the quickest.
#dumpGrid(vgrid, 'myfield.dat')       # binary, single precision
#dumpGridToTxt(vgrid, 'myfield.txt')  # ASCII

# load your own field
#loadGrid(vgrid, 'myfield.dat')
#loadGridToTxt(vgrid, 'myfield.txt')


# ### Running the simulation
# Now that we have our magnetic field ready we can fire up our simulation and hope that something visually interesting is going to happen.
# 
# Note that we use ParticleCollector to save the trajectory in memory. Alternatively, we could also use `TextOutput()`.

# In[3]:


sim = ModuleList()
sim.add(PropagationCK(Bfield))
sim.add(PhotoPionProduction(CMB))
sim.add(PhotoPionProduction(IRB))
sim.add(PhotoDisintegration(CMB))
sim.add(PhotoDisintegration(IRB))
sim.add(ElectronPairProduction(CMB))
sim.add(ElectronPairProduction(IRB))
sim.add(NuclearDecay())
sim.add(MaximumTrajectoryLength(25 * Mpc))
output = ParticleCollector()
output.setClone(True)
#output = TextOutput('trajectory.txt', Output.Trajectory3D) 
sim.add(output)

x = Vector3d(0,0,0)  # position
p = Vector3d(1,1,0)  # direction
c = Candidate(nucleusId(16, 8), 100 * EeV, x, p)

sim.run(c, True)


# ### (Optional) Plotting
# 
# We plot the trajectory of our oxygen-16 nucleus. To distinguish between secondary nuclei the following colors are used: protons are blue, alpha particles are green, everthing heavier is red.
# 
# If `TextOutput` is used, the trajectory data can be loaded with `data = genfromtxt('trajectory.txt', names=True)`.

# In[4]:


get_ipython().run_line_magic('matplotlib', 'inline')
from pylab import *
from mpl_toolkits.mplot3d import axes3d

data = {'X': [], 'Y': [], 'Z': [], 'ID': []}
for c in output:
    data['X'].append(c.current.getPosition().getX())
    data['Y'].append(c.current.getPosition().getY())
    data['Z'].append(c.current.getPosition().getZ())
    data['ID'].append(c.current.getId())

# trajectory points
x, y, z = np.array(data['X']), np.array(data['Y']), np.array(data['Z'])

# translate particle ID to charge number
Z = [chargeNumber(int(id)) for id in data['ID']]

# translate the charge number to color and size
# --> protons are blue, Helium is green, everthing else is red
colorDict = {0:'k', 1:'b', 2:'g', 3:'r', 4:'r', 5:'r', 6:'r', 7:'r', 8:'r'}
sizeDict = {0:4, 1:4, 2:8, 3:10, 4:10, 5:10, 6:10, 7:10, 8:10}
colors = [colorDict[z] for z in Z]
sizes  = [sizeDict[z] for z in Z]

fig = plt.figure(figsize=(14, 7))#plt.figaspect(0.5))
ax = fig.gca(projection='3d')# , aspect='equal'

ax.scatter(x/Mpc,y/Mpc,z/Mpc, 'o', s=sizes, lw=0, facecolor=colors)

ax.set_xlabel('x / Mpc', fontsize=18)
ax.set_ylabel('y / Mpc', fontsize=18)
ax.set_zlabel('z / Mpc', fontsize=18)
ax.set_xlim((-1, 16))
ax.set_ylim((-1, 16))
ax.set_zlim((-1, 16))
ax.xaxis.set_ticks((0, 5, 10, 15))
ax.yaxis.set_ticks((0, 5, 10, 15))
ax.zaxis.set_ticks((0, 5, 10, 15))

show()

