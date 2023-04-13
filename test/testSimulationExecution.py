import sys

try:
    import unittest
    import tempfile
except:
    print("***********************************************************")
    print("* WARNING!! Couldn't import python unittesting framework! *")
    print("* No python tests have been executed                      *")
    print("***********************************************************")
    sys.exit(0)

try:
    import crpropa as crp
except Exception as e:
    print("*** CRPropa import failed")
    print(type(e), str(e))
    sys.exit(-1)

class test1DChainWithSecondaries(unittest.TestCase):
    def runTest(self):

        outputFile = tempfile.NamedTemporaryFile()

        sim = crp.ModuleList()

        # photon fields
        CMB = crp.CMB()
        IRB = crp.IRB_Kneiske04()

        sim.add(crp.SimplePropagation(1 * crp.kpc, 1 * crp.Mpc))
        sim.add(crp.Redshift())
        sim.add(crp.PhotoPionProduction(CMB))
        sim.add(crp.PhotoPionProduction(IRB))
        sim.add(crp.PhotoDisintegration(CMB))
        sim.add(crp.PhotoDisintegration(IRB))
        sim.add(crp.NuclearDecay())
        sim.add(crp.ElectronPairProduction(CMB))
        sim.add(crp.ElectronPairProduction(IRB))
        sim.add(crp.MinimumEnergy(1 * crp.EeV))
        sim.add(crp.EMCascade())

        # observer
        obs = crp.Observer()
        obs.add(crp.Observer1D())
        sim.add(obs)

        # output
        output = crp.TextOutput(outputFile.name)
        output.set1D(True)
        obs.onDetection(output)

        # source
        source = crp.Source()
        source.add(crp.SourceUniform1D(1 * crp.Mpc, 1000 * crp.Mpc))
        source.add(crp.SourceRedshift1D())

        # power law spectrum with charge dependent maximum energy Z*100 EeV
        # elements: H, He, N, Fe with equal abundances at constant energy per
        # nucleon
        composition = crp.SourceComposition(1 * crp.EeV, 100 * crp.EeV, -1)
        composition.add(1,  1,  1)  # H
        composition.add(4,  2,  1)  # He-4
        composition.add(14, 7,  1)  # N-14
        composition.add(56, 26, 1)  # Fe-56
        source.add(composition)

        # run simulation
        sim.run(source, 1)

        outputFile.close()


if __name__ == '__main__':
    unittest.main()
