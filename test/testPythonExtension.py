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
    import crpropa as crp
except Exception as e:
    print("*** CRPropa import failed")
    print(type(e), str(e))
    sys.exit(-1)

import numpy as np



class testCrossLanguagePolymorphism(unittest.TestCase):

    def test_module(self):
        class CountingModule(crp.Module):
            def __init__(self):
                crp.Module.__init__(self)
                self.count = 0

            def process(self, c):
                self.count += 1

        count_accept = CountingModule()
        count_reject = CountingModule()
        filter = crp.ParticleFilter([-1, 1])
        filter.onAccept(count_accept)
        filter.onReject(count_reject)

        c = crp.Candidate()

        for id in [-1, 1, 6, 9, -19, 23, 100010001]:
            c.current.setId(id)
            filter.process(c)

    def test_ParticleCollector(self):
        c = crp.Candidate()
        p = crp.ParticleCollector()
        p.process(c)
        c_out = p[0]
        for c_i in p:
            c_out = c_i

    def test_ObserverFeature(self):
        class CountingFeature(crp.ObserverFeature):
            def __init__(self):
                crp.ObserverFeature.__init__(self)
                self.value = 0

            def checkDetection(self, candidate):
                self.value += 1
                return crp.DETECTED

        obs = crp.Observer()
        counter = CountingFeature()
        obs.add(counter)
        for i in range(5):
            candidate = crp.Candidate()
            obs.process(candidate)
            self.assertEqual(i + 1, counter.value)

    def testCustomMagneticField(self):
        class CustomMagneticField(crp.MagneticField):
            def __init__(self, val):
                crp.MagneticField.__init__(self)
                self.val = val
                
            def getField(self, position):
                return crp.Vector3d(self.val)

            def getField(self, position, z):
                return crp.Vector3d(self.val)

        field = CustomMagneticField(crp.gauss)
        propBP = crp.PropagationBP(field, 1e-4, 1*crp.Mpc, 1*crp.Mpc)
        propCK = crp.PropagationCK(field, 1e-4, 1*crp.Mpc, 1*crp.Mpc)
        propSDE = crp.DiffusionSDE(field)
        pos = crp.Vector3d(-1, 0, 0)
        z = 0
        fieldAtPos = field.getField(pos, z)
        self.assertEqual(fieldAtPos, propBP.getFieldAtPosition(pos, z))
        self.assertEqual(fieldAtPos, propCK.getFieldAtPosition(pos, z))
        self.assertEqual(fieldAtPos, propSDE.getMagneticFieldAtPosition(pos, z))

    def testCustomAdvectionField(self):
        class CustomAdvectionField(crp.AdvectionField):
            def __init__(self, val):
                crp.AdvectionField.__init__(self)
                self.val = val
                
            def getField(self, position):
                return crp.Vector3d(self.val)

            def getDivergence(self, position):
                return 0.0

        constMagVec = crp.Vector3d(0*crp.nG,0*crp.nG,1*crp.nG)
        magField = crp.UniformMagneticField(constMagVec)
        advField = CustomAdvectionField(1)
        propSDE = crp.DiffusionSDE(magField, advField)
        pos = crp.Vector3d(1, 0, 0)
        advFieldAtPos = advField.getField(pos)
        self.assertEqual(advFieldAtPos, propSDE.getAdvectionFieldAtPosition(pos))

    def testCustomMassDensity(self):
        class CustomMassDensity(crp.Density):
            def __init__(self, density, HIDensity, HIIDensity, H2Density, nucleonDensity):
                crp.Density.__init__(self)
                self.density = density
                self.HIDensity = HIDensity
                self.HIIDensity = HIIDensity
                self.H2Density = H2Density
                self.nucleonDensity = nucleonDensity
                
            def getDensity(self, position):
                return self.density

            def getHIDensity(self, position):
                return self.HIDensity

            def getHIIDensity(self, position):
                return self.HIIDensity

            def getH2Density(self, position):
                return self.H2Density

            def NucleonDensity(self, position):
                return self.nucleonDensity

        density = 10
        HIDensity = 5
        HIIDensity = 2
        H2Density = 1
        nucleonDensity = 0.5
        massDensity = CustomMassDensity(density, HIDensity, HIIDensity, H2Density, nucleonDensity)
        pos = crp.Vector3d(1, 0, 0)
        self.assertEqual(density, massDensity.getDensity(pos))
        self.assertEqual(HIDensity, massDensity.getHIDensity(pos))
        self.assertEqual(HIIDensity, massDensity.getHIIDensity(pos))
        self.assertEqual(H2Density, massDensity.getH2Density(pos))

    def testCustomPhotonField(self):
        class CustomPhotonField(crp.PhotonField):
            def __init__(self, density):
                crp.PhotonField.__init__(self)
                self.density = density
                self.fieldName = 'testCustomPhotonField'
                self.isRedshiftDependent = True
                
            def getPhotonDensity(self, energy, z):
                return self.density

            def getFieldName(self):
                return self.fieldName

            def hasRedshiftDependence(self):
                return self.isRedshiftDependent

        photonDensity = 10
        photonField = CustomPhotonField(photonDensity)
        energy = 10*crp.GeV
        z = 0
        self.assertEqual(photonDensity, photonField.getPhotonDensity(energy, z))
        self.assertEqual('testCustomPhotonField', photonField.getFieldName())
        self.assertEqual(True, photonField.hasRedshiftDependence())


    def testCustomPhotonField(self):
        class CustomPhotonField(crp.PhotonField):
            def __init__(self, val):
                crp.PhotonField.__init__(self)
                self.val = val

            def getFieldName(self):
                return 'CMB'

            def getPhotonDensity(self, ePhoton, z):
                return self.val
                
        constDensity = 1
        photonField = CustomPhotonField(constDensity)
        ppp = crp.PhotoPionProduction(photonField)
        pppPhotonField = ppp.getPhotonField()

        self.assertEqual(constDensity, pppPhotonField.getPhotonDensity(0))


class testCandidatePropertymap(unittest.TestCase):
    def setUp(self):
        self.candidate = crp.Candidate()

    def __propertySetGet(self, value):
        self.candidate.setProperty('Foo', value)
        self.assertEqual(value, self.candidate.getProperty('Foo'))

    def testString(self):
        self.__propertySetGet('Bar')

    def testUnicode(self):
        self.__propertySetGet(u'Bar')

    def testUnicodeName(self):
        self.candidate.setProperty(u'Foo', 23)
        self.assertEqual(23, self.candidate.getProperty(u'Foo'))

    def testBool(self):
        self.__propertySetGet(True)
        self.__propertySetGet(False)

    def testInt(self):
        self.__propertySetGet(42)

    def testFloat(self):
        self.__propertySetGet(3.14)
        v = np.array([2.])
        self.__propertySetGet(v[0])


class testKeywordArguments(unittest.TestCase):
  def testExceptionOnNonExistingArguemnt(self):
    with self.assertRaises(Exception, msg="This is likely due to a swig bug. Please try to disable the builtin option by compiling crpropa with cmake .. -DENABLE_SWIG_BUILTIN=OFF"):
      p = crp.PhotoDisintegration(nonExistingKeywordArguemntShouldRaiseException=True)
  def testDisablingOfKwargs(self):
    with self.assertRaises(Exception, msg="This is likely due to a swig bug. Please try to disable the builtin option by compiling crpropa with cmake .. -DENABLE_SWIG_BUILTIN=OFF"):
      p = crp.PhotoDisintegration(photonField=crp.IRB_Dominguez11)
  # swig currently does not support kwargs in overloaded functions - we should
  # thus disable them.
  #def testKeywordArgument(self):
  #  p = crp.PhotoDisintegration(photonField=crp.IRB_Dominguez11)
  #  self.assertTrue('IRB_Dominguez11' in p.getDescription())


class testVector3(unittest.TestCase):
  def testPublicReferenceAccess(self):
    v = crp.Vector3d(1., 2., 3.)
    self.assertEqual(v.x, 1.)
    self.assertEqual(v.y, 2.)
    self.assertEqual(v.z, 3.)
    v.x = 23.
    self.assertEqual(v.x, 23.)

    ## this test fails in some systems
    # def testArrayInterface(self):
    #   # this test fails for some combinations of Python version and system
    #   v = crp.Vector3d(1., 2., 3.)
    #   self.assertEqual(2., np.mean(v) )
    #   x = np.ones(3)
    #   self.assertEqual(6., sum(v * x) )

  def testRepr(self):
    v = crp.Vector3d(1., 2., 3.)
    if sys.version_info >= (3, 3):
      import unittest.mock
      import io
      with unittest.mock.patch('sys.stdout', new = io.StringIO()) as fake_out:
        print(v)
    else:
      import StringIO
      fake_out = StringIO.StringIO()
      sys.stdout = fake_out
      print(v)
      sys.stdout = sys.__stdout__
    self.assertEqual(fake_out.getvalue().rstrip(), v.getDescription())

  def testOutOfBound(self):
    v = crp.Vector3d(1., 2., 3.)
    self.assertRaises(IndexError, v.__getitem__, 3)
    self.assertRaises(IndexError, v.__setitem__, 3, 10)

class testParticleCollector(unittest.TestCase):
  def testParticleCollectorIterator(self):
    collector = crp.ParticleCollector()
    lengths = [1*crp.pc, 10*crp.pc, 100*crp.pc]
    for l in lengths:
        c = crp.Candidate()
        c.setTrajectoryLength(l)
        collector.process(c)

    self.assertEqual(len(collector), len(lengths))

    for c, l in zip(collector, lengths):
        self.assertEqual(c.getTrajectoryLength(), l)

  def testParticleCollectorAsModuleListInput(self):
    sim = crp.ModuleList()
    sim.add(crp.MaximumTrajectoryLength(3.14))
    sim.add(crp.SimplePropagation(0.001, 0.001))
    collector = crp.ParticleCollector()
    c1 = crp.Candidate()
    c2 = crp.Candidate()
    collector.process(c1)
    collector.process(c2)
    sim.run(collector.getContainer())
    for c in collector:
        self.assertAlmostEqual(
            c.getTrajectoryLength(), 3.14, places=2)

  def testParticleCollectorAsModuleListOutput(self):
    sim = crp.ModuleList()
    sim.add(crp.MaximumTrajectoryLength(3.14))
    sim.add(crp.SimplePropagation(0.001, 0.001))
    collector = crp.ParticleCollector()
    sim.add(collector)
    c = crp.Candidate()
    sim.run(c)
    self.assertAlmostEqual(
        collector[0].getTrajectoryLength(),
        3.14, places=2)

class testGrid(unittest.TestCase):
  def testGridPropertiesConstructor(self):
    N = 32
    gp = crp.GridProperties(crp.Vector3d(0), N, 0.1)
    grid = crp.Grid1f(gp)
    self.assertEqual(grid.getNx(), 32)

if hasattr(crp, 'GridTurbulence'):
    class testTurbulentField(unittest.TestCase):
      #check problems brought up in https://github.com/CRPropa/CRPropa3/issues/322
      def testTurbulenceSpectrum(self):
        spectrum = crp.TurbulenceSpectrum(1., 1., 10.)
        self.assertEqual(spectrum.getBrms(), 1.)
        self.assertEqual(spectrum.getLmin(), 1.)
        self.assertEqual(spectrum.getLmax(), 10.)
        self.assertEqual(spectrum.getLbendover(), 1.)
        self.assertEqual(spectrum.getSindex(), 5./3.)
        self.assertEqual(spectrum.getQindex(), 4.)

      def testGridTurbulence(self):
        N = 64
        boxSize = 1*crp.Mpc
        l_bo = boxSize/8
        spacing = boxSize / N
        tf = crp.GridTurbulence(
            crp.TurbulenceSpectrum(1.0, 2*spacing, boxSize, l_bo),
            crp.GridProperties(crp.Vector3d(0), N, spacing)
        )

if __name__ == '__main__':
    unittest.main()
