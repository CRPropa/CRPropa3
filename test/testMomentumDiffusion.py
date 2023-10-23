# coding=utf-8
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
    print("* WARNING!! Couldn't import numpy framework! *")
    print("* No python tests have been executed                      *")
    print("***********************************************************")
    sys.exit(-1)

try:
    import crpropa
    from crpropa import nG, kpc, pc, GeV, TeV, PeV, c_light
except Exception as e:
    print("*** CRPropa import failed")
    print(type(e), str(e))
    sys.exit(-2)


class MomentumDiffusion(unittest.TestCase):
	
    Dpp = 10
    limit = 0.3
	
    momDif = crpropa.ConstantMomentumDiffusion(Dpp)
    momDif.setLimit(limit)
	
    def test_Simple(self):
        self.assertEqual(self.momDif.getDpp(), self.Dpp)
        self.assertEqual(self.momDif.getLimit(), self.limit)

    def test_Limits(self):
        Dpp2 = 10
        limit2 = 0.4
	
        momDif = crpropa.ConstantMomentumDiffusion(Dpp2, limit2)
        self.assertEqual(momDif.getDpp(), Dpp2)
        self.assertEqual(momDif.getLimit(), limit2)
    
    def test_Helper(self):
        """Test to check the calculation of the helper functions
    	Dpp
    	B
    	A
    	"""
        E = 10*TeV
        c = crpropa.Candidate(crpropa.nucleusId(1,1))
        c.current.setEnergy(E)
        p = c.current.getEnergy() / c_light
        
        bscal = self.momDif.calculateBScalar()
        self.assertEqual(np.sqrt(2*self.Dpp), bscal)
        ascal = self.momDif.calculateAScalar(p, Dpp)
        self.assertAlmostEqual(2/p*self.Dpp, ascal)
    
    def test_NeutralParticle(self):
        E = 10*TeV
        c = crpropa.Candidate(crpropa.nucleusId(1,0))
        c.current.setEnergy(E)
        c.setNextStep(10)
        self.momDif.process(c)
        
        self.assertEqual(c.current.getEnergy(), E) #acts only on charged particles
        self.assertEqual(c.getNextStep(), 10)
