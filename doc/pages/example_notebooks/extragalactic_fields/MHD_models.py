
# coding: utf-8

# # 3D MHD models
# This notebook explains how to use cubic results of 3D MHD models on a uniform grid in CRPropa.
# 
# ## Supplied data
# 
# The fields need to be supplied in a raw binary file that contains only single floats, arranged as follows: Starting with the cell values (Bx,By,Bz for magnetic field or rho for density) at the origin of the box, the code continues to read along z, then y and finally x.
# 
# On https://crpropa.desy.de/ under "Additional resources" you can find a number of MHD models used with CRPropa in the literature. 
# 
# 

# In[1]:

from crpropa import *

## settings for MHD model (must be set according to model)
filename_bfield = "clues_primordial.dat" ## filename of the magnetic field
gridOrigin = Vector3d(0,0,0)             ## origin of the 3D data, preferably at boxOrigin
gridSize = 1024                          ## size of uniform grid in data points
size = 249.827*Mpc                       ## physical edgelength of volume in Mpc
b_factor = 1.                            ## global renormalizatino factor for the field

## settings of simulation
boxOrigin = Vector3d( 0, 0, 0,)          ## origin of the full box of the simulation
boxSize = Vector3d( size, size, size )   ## end of the full box of the simulation

## settings for computation
minStep = 10.*kpc                        ## minimum length of single step of calculation
maxStep = 4.*Mpc                         ## maximum length of single step of calculation
tolerance = 1e-2                         ## tolerance for error in iterative calculation of propagation step

spacing = size/(gridSize)                ## resolution, physical size of single cell

m = ModuleList()


## instead of  computing propagation without Lorentz deflection via
# m.add(SimplePropagation(minStep,maxStep))

## initiate grid to hold field values
vgrid = VectorGrid( gridOrigin, gridSize, spacing )
## load values to the grid
loadGrid( vgrid, filename_bfield, b_factor )
## use grid as magnetic field
bField = MagneticFieldGrid( vgrid )
## add propagation module to the simulation to activate deflection in supplied field
m.add(PropagationCK( bField, tolerance, minStep, maxStep))
#m.add(DeflectionCK( bField, tolerance, minStep, maxStep))  ## this was used in older versions of CRPropa


# to make use of periodicity of the provided data grid, use

# In[2]:

m.add( PeriodicBox( boxOrigin, boxSize ) )


# to not follow particles forever, use

# In[3]:

m.add( MaximumTrajectoryLength( 400*Mpc ) ) 


# ## Uniform injection
# 
# The most simple scenario of UHECR sources is a uniform distribution of their sources. This can be realized via use of

# In[4]:

source = Source()
source.add( SourceUniformBox( boxOrigin, boxSize )) 


# ## Injection following density field
# 
# The distribution of gas density can be used as a probability density function for the injection of particles from random positions.

# In[65]:

filename_density = "mass-density_clues.dat" ## filename of the density field

source = Source()
## initialize grid to hold field values
mgrid = ScalarGrid( gridOrigin, gridSize, spacing )
## load values  to grid
loadGrid( mgrid, filename_density )
## add source module to simulation
source.add( SourceDensityGrid( mgrid ) )


# ## Mass Halo injection
# 
# Alternatively, for the CLUES models, we also provide a list of mass halo positions. These positions can be used as sources with same properties by use of the following

# In[67]:

filename_halos = 'clues_halos.dat'

# read data from file
data = np.loadtxt(filename_halos, unpack=True, skiprows=39)
sX = data[0]                                                                                 
sY = data[1]                                                                                 
sZ = data[2]                                                                                 
mass_halo = data[5]                                                                          

## find only those mass halos inside the provided volume (see Hackstein et al. 2018 for more details)
Xdown= sX >= 0.25                                                                            
Xup= sX <= 0.75                                                                              
Ydown= sY >= 0.25                                                                            
Yup= sY <= 0.75                                                                              
Zdown= sZ >= 0.25                                                                            
Zup= sZ <= 0.75                                                                              
insider= Xdown*Xup*Ydown*Yup*Zdown*Zup                                                       

## transform relative positions to physical positions within given grid
sX = (sX[insider]-0.25)*2*size
sY = (sY[insider]-0.25)*2*size
sZ = (sZ[insider]-0.25)*2*size

## collect all sources in the multiple sources container
smp = SourceMultiplePositions()
for i in range(0,len(sX)):
    pos = Vector3d( sX[i], sY[i], sZ[i] )
    smp.add( pos, 1. )
    
## add collected sources
source = Source()
source.add( smp )


# additional source properties

# In[5]:

## use isotropic emission from all sources
source.add( SourceIsotropicEmission() )

## set particle type to be injected
A, Z = 1, 1 # proton
source.add( SourceParticleType( nucleusId(A,Z) ) )

## set injected energy spectrum
Emin, Emax = 1*EeV, 1000*EeV
specIndex = -1
source.add( SourcePowerLawSpectrum( Emin, Emax, specIndex ) ) 


# ## Observer
# 
# To register particles, an observer has to be defined. In the provided constrained simulations the position of the Milky Way is, by definition, in the center of the volume.

# In[6]:

filename_output = 'output.txt'

filename_output = 'data/output_MW.txt'

obsPosition = Vector3d(0.5*size,0.5*size,0.5*size) # position of observer, MW is in center of constrained simulations
obsSize = 800*kpc  ## physical size of observer sphere


## initialize observer that registers particles that enter into sphere of given size around its position
obs = Observer()
obs.add( ObserverSmallSphere( obsPosition, obsSize ) )
## write registered particles to output file
obs.onDetection( TextOutput( filename_output ) )
## choose to not further follow particles paths once detected
obs.setDeactivateOnDetection(True)
## add observer to module list
m.add(obs)


# finally run the simulation by

# In[7]:

N = 1000

m.showModules()         ## optional, see summary of loaded modules
m.setShowProgress(True) ## optional, see progress during runtime
m.run(source, N, True)  ## perform simulation with N particles injected from source

