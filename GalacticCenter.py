from crpropa import *
from pylab import *
import random as rd
import numpy as np
import math
HI=HadronicInteraction()
#Number of particles
particles = 10000


#Candidate
c = Candidate()

#Radius of Sphere

r= 234*pc

#Steplength
steplength = 0.001

#Initialise Sourcelist
sl=SourceList()
origin = Vector3d(0)

#Create Sources
#Calculated by Integeral
for i in range (50):
    s=Source()
    s.add(SourcePosition(HI.Position(60*pc,200*pc)))
    s.add(SourceParticleType(nucleusId(1,1)))
    s.add(SourcePowerLawSpectrum(1*TeV, 10000*TeV, -2.0))
    s.add(SourceIsotropicEmission())
    sl.add(s)

#Create output1
output = TextOutput("../../Simulations310518/1_ObS.txt")
output.disableAll()
output.enable(output.SerialNumberColumn)
output.enable(output.CurrentIdColumn)
output.enable(output.CurrentEnergyColumn)
output.enable(output.CurrentPositionColumn)
output.enable(output.SourcePositionColumn)
output.setEnergyScale(TeV)
output.setLengthScale(pc)

#Create Observer Sphere
obs=Observer()
obs.add(ObserverLargeSphere(origin, r))
#obs.setDeactivateOnDetection(False)
obs.onDetection(output)


#Create output2
output2 = TextOutput("../../Simulations310518/1_Source.txt")
output2.disableAll()
output2.enable(output2.SerialNumberColumn)
output2.enable(output2.CurrentIdColumn)
output2.enable(output2.CurrentEnergyColumn)
output2.enable(output2.CurrentPositionColumn)
output2.enable(output2.SourcePositionColumn)
output2.enable(output2.TrajectoryLengthColumn)
output2.setEnergyScale(TeV)
output2.setLengthScale(pc)

#Create Time Evolution Source
obs2=Observer()
obs2.add(ObserverTimeEvolution(0, 1* pc, 1))
obs2.setDeactivateOnDetection(False)
obs2.onDetection(output2)

#Create output4  kpc
output4 = TextOutput("../../Simulations310518/1_Traj3kpc.txt")
output4.disableAll()
output4.enable(output4.SerialNumberColumn)
output4.enable(output4.CurrentIdColumn)
output4.enable(output4.CurrentEnergyColumn)
output4.enable(output4.CurrentPositionColumn)
output4.enable(output4.SourcePositionColumn)
output4.enable(output4.TrajectoryLengthColumn)
output4.setEnergyScale(TeV)
output4.setLengthScale(pc)

#Create Time Evolution kpc
obs4=Observer()
obs4.add(ObserverTimeEvolution(0*kpc, 0*kpc, 2))
obs4.setDeactivateOnDetection(False)
obs4.onDetection(output4)






#Initiate turbulent magnetic field
MagField = MagneticFieldList()
randomSeed = 23
lMin=0.04 * pc
lMax=2.*pc
l= turbulentCorrelationLength(lMin, lMax, -11./3.)
spacing=0.01*pc
vgrid = VectorGrid(origin, 600, spacing)
vgrid.setReflective(True)
b=50*1e-6*gauss
initTurbulence(vgrid, b, lMin, lMax, -11./3., randomSeed)
TurbField = MagneticFieldGrid(vgrid)
MagField.addField(TurbField)

#GC=GalacticCenter()


#vgrid.setReflective(True)
#JF.setTurbulentGrid(vgrid)
#JF.setStriatedGrid(vgrid)
#JF.setRegularGrid(vgrid)

MagField.addField(TurbField)

#ModuleList
m = ModuleList()
#Propagation
mod_Diffision=DiffusionSDE(MagField, 0.5, steplength * pc, 100 * steplength * pc)
m.add(mod_Diffision)

#Interaction
m.add(HadronicInteraction())
m.add(EMInverseComptonScattering(Blackbody, True, 1))
m.add(EMDoublePairProduction(Blackbody, True, 1))
m.add(EMPairProduction(Blackbody, True, 1))
m.add(EMTripletPairProduction(Blackbody, True, 1))
m.add(ElasticScattering(Blackbody))
m.add(ElectronPairProduction(Blackbody, True, 1))
m.add(SynchrotronRadiation(TurbField, True, 1))
m.add(PhotoDisintegration(Blackbody, True, 1))
m.add(PhotoPionProduction(Blackbody, True, True, False, 1, False))

#Observer
m.add(obs)
m.add(obs2)

m.add(obs4)

m.add(MaximumTrajectoryLength(5*kpc))

m.setShowProgress(True)
m.run(sl, particles, True)

#Explicitly closing files.
#Important e.g. for the usage with jupyter notebooks
output.close()
output2.close()
output4.close()






