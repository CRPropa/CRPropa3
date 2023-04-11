import crpropa
import myPlugin

print("My Simulation\n")

ml = crpropa.ModuleList()

ml.add(crpropa.SimplePropagation(1*crpropa.parsec, 100*crpropa.parsec))
ml.add(crpropa.MaximumTrajectoryLength(1000*crpropa.parsec))
ml.add(myPlugin.MyModule())

print("+++ List of modules")
print(ml.getDescription())


print("+++ Preparing source")
source = crpropa.Source()
source.add(myPlugin.AddMyProperty())
print(source.getDescription())

print("+++ Starting Simulation")
ml.run(source, 1)

print("+++ Done")
