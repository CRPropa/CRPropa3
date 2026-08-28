"""
Python-side tests for PhotoPionProductionEmpirical.

The physics of the empirical photomeson model (Morejon et al., JCAP 11 (2019)
007) is pinned in C++, in testInteraction.cpp.  This file covers two elements
not covered by the C++ tests:

  * The SWIG bindings: check that the class constructs from Python, is visible
    as a PhotoPionProduction / Module subclass, and exposes its elements.

  * The module called inside an example propagation loop. check it runs the 
    full chain of rate lookup -> distance sampling -> eps sampling ->
    interaction -> secondary collection, repeated over a trajectory.

TestPhysics deliberately repeats a few assertions from the C++ suite, but just
to verify the exposed functions work correctly through the bindings.

The data tables are produced by the CRPropa3-data repository. When they are 
missing the module raises a RuntimeError pointing at 
doc/photomesonEmpirical/README.md, which documents where they come from and
how to point CRPROPA_DATA_PATH at them; every test here turns that into a skip,
so a checkout without the tables reports skipped rather than failed.
"""

import sys
import unittest

try:
    import crpropa as crp
except Exception as e:
    print("*** CRPropa import failed:", e)
    sys.exit(-1)

CMB = crp.CMB()


def _module(*args):
    try:
        return crp.PhotoPionProductionEmpirical(CMB, *args)
    except RuntimeError as e:
        raise unittest.SkipTest("empirical photomeson tables not generated: %s" % e)


class TestConstruction(unittest.TestCase):
    def test_is_a_photopion_module(self):
        m = _module()
        self.assertIsInstance(m, crp.PhotoPionProduction)
        self.assertIsInstance(m, crp.Module)
        self.assertEqual(m.getInteractionTag(), "PPPE")

    def test_settings_round_trip(self):
        m = _module()
        m.setDecayMode(crp.PhotoPionProductionEmpirical.DisplacedDecay)
        self.assertEqual(m.getDecayMode(),
                         crp.PhotoPionProductionEmpirical.DisplacedDecay)
        m.setMultiplicitySampling(crp.PhotoPionProductionEmpirical.PoissonSampling)
        self.assertEqual(m.getMultiplicitySampling(),
                         crp.PhotoPionProductionEmpirical.PoissonSampling)
        m.setHaveKaons(False)
        self.assertFalse(m.getHaveKaons())
        m.setLimit(0.05)
        self.assertAlmostEqual(m.getLimit(), 0.05)

    def test_rate_needs_no_tabulated_photon_field(self):
        # the rate is built from the field itself, so any field not
        # included in CRPropa's tabulated photopion rates still work
        m = _module()
        m.setPhotonField(crp.BlackbodyPhotonField("custom", 10.0))
        self.assertGreater(m.getRate(1, 1, 1e12, 0), 0)


class TestPhysics(unittest.TestCase):
    def test_proton_rate_matches_sophia(self):
        # for A = 1 the empirical model matches SOPHIA
        m = _module()
        ppp = crp.PhotoPionProduction(CMB)
        for lg in (11.0, 11.5, 12.0, 12.5, 13.0):
            gamma = 10 ** lg
            self.assertAlmostEqual(
                m.getRate(1, 1, gamma, 0) * ppp.nucleonMFP(gamma, 0, True), 1.0,
                delta=0.05, msg="log10(gamma) = %.1f" % lg)

    def test_iron_shadowing(self):
        # above 1.2 GeV: sigma ~ A^alpha with alpha -> 0.66
        m = _module()
        eps = m.getEnergyBin(140)
        superposition = 26 * m.crossection(eps, 1, 1) + 30 * m.crossection(eps, 1, 0)
        self.assertAlmostEqual(m.crossection(eps, 56, 26) / superposition,
                               56 ** -0.34, delta=1e-3)

    def test_iron_pion_scaling(self):
        m = _module()
        self.assertLess(m.pionScaling(m.getEnergyBin(20), 56, 26), 0.6)
        self.assertAlmostEqual(m.pionScaling(m.getEnergyBin(140), 56, 26), 1.0,
                               delta=1e-3)


class TestSimulation(unittest.TestCase):
    # Photopion production needs log10(gamma) >~ 11 on the CMB, and for a nucleus
    # gamma ~ E / (A m_N): Fe-56 has to carry several thousand EeV before it
    # interacts at all.  Unphysical as a source spectrum, but acceptable for a 
    # test of the module.
    ENERGY = 20000 * crp.EeV

    def _run(self, module, source_id, n=200):
        sim = crp.ModuleList()
        sim.add(crp.SimplePropagation(1 * crp.kpc, 10 * crp.Mpc))
        sim.add(module)
        sim.add(crp.MinimumEnergy(1 * crp.EeV))
        sim.add(crp.MaximumTrajectoryLength(500 * crp.Mpc))
        source = crp.Source()
        source.add(crp.SourceParticleType(source_id))
        source.add(crp.SourceEnergy(self.ENERGY))
        source.add(crp.SourcePosition(crp.Vector3d(0)))
        source.add(crp.SourceDirection(crp.Vector3d(1, 0, 0)))
        collector = crp.ParticleCollector()
        sim.add(collector)
        sim.setShowProgress(False)
        sim.run(source, n, True)
        return collector

    def test_iron_breaks_up_much_more_than_superposition(self):
        iron = crp.nucleusId(56, 26)

        def mean_surviving_mass(module):
            collector = self._run(module, iron)
            masses = [crp.massNumber(c.current.getId())
                      for c in collector
                      if crp.isNucleus(c.current.getId())
                      and crp.massNumber(c.current.getId()) > 4]
            return sum(masses) / float(len(masses)) if masses else 56.0

        empirical = mean_surviving_mass(_module(True, True, True))
        sophia = mean_surviving_mass(crp.PhotoPionProduction(CMB, True, True, True))
        # both are photodisintegration-free, so any difference here 
        # comes from the photomeson mass loss in the empirical model,
        # which is much stronger than in SOPHIA
        self.assertLess(empirical, sophia)

    def test_neutrinos_are_produced(self):
        collector = self._run(_module(True, True, True), crp.nucleusId(56, 26))
        neutrinos = [c for c in collector
                     if abs(c.current.getId()) in (12, 14)]
        self.assertGreater(len(neutrinos), 0)


if __name__ == "__main__":
    unittest.main()
