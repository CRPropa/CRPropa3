from crpropa import *
from pylab import *
import random as rd
import numpy as np
import math
import pandas as pd

#Number of particles
particles = 10

#Radius of Sphere
r=162*pc

#Steplength
steplength = 0.1

#Initialise Sourcelist
sl=SourceList()
origin = Vector3d(0)

#Create 2600 Sources
for i in range (2): # Attention: source number reduced by a factor of 100
    s=Source()
    phi = 2*M_PI*rd.random()
    RandRadius = 150*pc*pow(rd.random(), 1. / 2.)
    a= Vector3d (np.cos(phi)*RandRadius, np.sin(phi)*RandRadius, (-0.5+rd.random())*120*pc)
    s.add(SourcePosition(a))
    s.add(SourceParticleType(nucleusId(1,0)))
    #s.add(SourceParticleType(14))
    s.add(SourcePowerLawSpectrum(1*TeV, 10*TeV, - 2.0)) # hier sollte "minus" 2 stehen
    s.add(SourceIsotropicEmission()) # sonst fliegen die alle in die gleiche Richtung los
    sl.add(s)

#Create output1
output = TextOutput("OS.txt")
output.disableAll()
output.enable(Output.SerialNumberColumn)
output.enable(Output.CurrentIdColumn)
output.enable(Output.CurrentEnergyColumn)
output.enable(Output.CurrentPositionColumn)
output.enable(Output.SourcePositionColumn)
output.setEnergyScale(TeV)
output.setLengthScale(pc)

#Create Observer Sphere
obs=Observer()
obs.add(ObserverLargeSphere(Vector3d(0.), r))
obs.onDetection(output)

#Create output2
output2 = TextOutput("SourceTO.txt")
output2.disableAll()
output2.enable(Output.SerialNumberColumn)
output2.enable(Output.CurrentIdColumn)
output2.enable(Output.CurrentEnergyColumn)
output2.enable(Output.CurrentPositionColumn)
output2.enable(Output.SourcePositionColumn)
output2.setEnergyScale(TeV)
output2.setLengthScale(pc)

#Create Time Evolution (original and 5 kpc)
obs2=Observer()
obs2.add(ObserverTimeEvolution(1*pc, 5*kpc, 5))
obs2.setDeactivateOnDetection(False)
obs2.onDetection(output2)

#Initiate turbulent magnetic field
MagField = MagneticFieldList()
randomSeed = 23
lMin=0.04 * pc
lMax=2.*pc
l= turbulentCorrelationLength(lMin, lMax, -11./3.)
spacing=0.01*pc
vgrid = VectorGrid(origin, 600, spacing)
vgrid.setReflective(True)
b=350*1e-6*gauss
initTurbulence(vgrid, b, lMin, lMax, -11./3., randomSeed)
TurbField = MagneticFieldGrid(vgrid)
MagField.addField(TurbField)

#ModuleList
m = ModuleList()
#Propagation
mod_PropCK=PropagationCK(MagField, 1e-3, steplength * pc, 100 * steplength * pc)
m.add(mod_PropCK)

#Interaction
m.add(HadronicInteraction())
m.add(EMInverseComptonScattering(Blackbody, True, True))
m.add(EMDoublePairProduction(Blackbody, True, True))
m.add(EMPairProduction(Blackbody, True, True))
m.add(EMTripletPairProduction(Blackbody, True, True))
m.add(ElasticScattering(Blackbody))
m.add(ElectronPairProduction(Blackbody, True, True))
m.add(SynchrotronRadiation(TurbField, True, True))
m.add(PhotoDisintegration(Blackbody, True, True))
m.add(PhotoPionProduction(Blackbody, True, True, False, True, False))

#Observer
m.add(obs)
m.add(obs2)
m.add(MaximumTrajectoryLength(10*kpc))

m.setShowProgress(True)
m.run(sl, particles, True)

#Explicitly closing files.
#Important e.g. for the usage with jupyter notebooks
output.close()
output2.close()


