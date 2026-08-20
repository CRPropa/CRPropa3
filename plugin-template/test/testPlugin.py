import sys

try:
    import unittest
except:
    print("***********************************************************")
    print("* WARNING!! Couldn't import python unittesting framework! *")
    print("* No python tests have been executed                      *")
    print("***********************************************************")
    sys.exit(0)

try:
    import numpy as np
except:
    print("***********************************************************")
    print("* WARNING!! Couldn't import numpy framework!              *")
    print("* No python tests have been executed                      *")
    print("***********************************************************")
    sys.exit(-1)

try:
    import crpropa
except Exception as e:
    print("*** CRPropa import failed")
    print(type(e), str(e))
    sys.exit(-2)
    
try:
    import myPlugin
except Exception as e:
    print("*** myPlugin import failed")
    print(type(e), str(e))
    sys.exit(-2)

class myPluginTests(unittest.TestCase):
	def testSimpleSimulation(self):
		ml = crpropa.ModuleList()

		ml.add(crpropa.SimplePropagation(1 * crpropa.parsec, 100 * crpropa.parsec))
		ml.add(crpropa.MaximumTrajectoryLength(1000 * crpropa.parsec))
		ml.add(myPlugin.MyModule())

		# Preparing source
		source = crpropa.Source()
		source.add(myPlugin.AddMyProperty())

		# Starting Simulation
		ml.run(source, 1)
		# Done


if __name__ == '__main__':
    unittest.main()