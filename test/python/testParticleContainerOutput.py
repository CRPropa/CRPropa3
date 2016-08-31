from crpropa import *

sim = ModuleList()
sim.add(SimplePropagation(1*kpc, 10*Mpc))

obs = Observer()
obs.add(ObserverPoint())
output = ParticleContainerOutput()
obs.onDetection(output)
sim.add(obs)

source = Source()
source.add(SourceUniform1D(1 * Mpc, 1000 * Mpc))
source.add(SourceParticleType(nucleusId(1,1)))

# run simulation
sim.setShowProgress(True)
sim.run(source, 20000, True)

for c in output:
    print c
