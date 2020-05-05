The CRPropa functionality can be extended in several ways using the provided
python interface.

### Implementing additional Modules and Features using Python only
Derivatives of base classes such as Module, Source, SourceProperty, etc. can be implemented in
Python and used like the built-in classes. Here is an example for a custom
Module:
```python
class MyModule(Module):
    """ Reduces the cosmic ray energy by 10% in each step """
    def process(self, c):
        c.current.setEnergy(c.current.getEnergy() * 0.9)

m = ModuleList()

mod = MyModule() # See https://github.com/CRPropa/CRPropa3/issues/165
m.add(mod)

c = Candidate()
c.current.setEnergy(10)
m.process(c)
print c.current.getEnergy()
```

When redefining the constructor make sure to call the constructor of the super
classes as well, as otherwise the code will segfault.
```python
class MyModule(Module):
    def __init__(self):
        Module.__init__(self)
```


The initial properties of a cosmic rays can be set with a Source, composed of several SourceProperties.
Custom SourceFeatures can be written in the following way:
```python
class MySourceFeature(SourceFeature):
    """ Set the initial energy to 10 EeV """
    def __init__(self):
        SourceFeature.__init__(self)

    def prepareParticle(self, particleState):
        particleState.setEnergy(10 * EeV)

s = Source()
s.add(MySourceFeature())
c = s.getCandidate()
print c.current.getEnergy()
```

The redshift is stored in the Candidate, not in the ParticleState. To set it with a SourceFeature use the following:
```python
class MySourceFeature(SourceFeature):
    """ Set the initial redshift """
    def __init__(self):
        SourceFeature.__init__(self)

    def prepareCandidate(self, candidate):
        candidate.setRedshift(0.6)

# The source feature has to be created outside of the class attribute
# s.add(MySourceFeature()) wil NOT work!
srcFtr = MySourceFeature()
s = Source()
s.add(srcFtr)
c = s.getCandidate()
print c.getRedshift()
```


### Manual Simulation Processing
If necessary, the simulation chain (ModuleList) and the cosmic ray source (Source, SourceFeature) can be
replaced by custom loops.

```python
m1 = SimplePropagation()
m2 = ElectronPairProduction()
m3 = MinimumEnergy(5 * EeV)

while mycandidate.isActive():
    m1.process(mycandidate)
    m2.process(mycandidate)
    m3.process(mycandidate)
    # reduce the energy by 10%
    E = mycandidate.current.getEnergy()
    mycandidate.current.setEnergy(0.9 * E)
```

A source can be replaced by setting all necessary cosmic rays properties by hand.
The created Candidates can either be propagated one-by-one, or first collected
in a CandidateVector and then propagated.
```python
for i in range(1000):
    p = ParticleState()
    p.setId(nucleusId(12, 6))
    p.setEnergy(200 * EeV)
    p.setPosition(Vector3d(100, 10, 10) * Mpc)
    p.setDirection(Vector3d(-1,0,0))

    c = Candidate(p)
    m.process(c)
```


### Plugins: Integrate Custom C++ Code to CRPropa's Python Steering
Extending CRPropa with C++ code and keep python steering is also possible using
SWIG.  This allows to integrate your code seamless as e.g.
```
import crpropa
import myPlugin

ml = crpropa.ModuleList()
ml.add(crpropa.MaximumTrajectoryLength(1000*crpropa.parsec))
ml.add(myPlugin.MyModule())

source = crpropa.Source()
source.add(myPlugin.AddMyProperty())
```
A template  is in the [plugin-template
folder](https://github.com/CRPropa/CRPropa3/tree/master/plugin-template) of the
CRPropa source. Although the template is complete with build and SWIG wrapper
code, deeper knowledge of SWIG, C++, and CMake are likely required for complex
projects.

To get started with your own plugin

1. Copy the folder to a new location. We highly recommended to manage the files
in a (git) repository from the beginning.
2. Test compiling the template
```
mkdir build
cd build
cmake ..
make && python ../testPlugin.py
```
This should work if CRPropa is installed and can be found by python.

3. Customize the template to your needs, starting with
naming your plugin by modifying the according line in `CMakeLists.txt' and
renaming the files accordingly.
