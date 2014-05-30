CRPropa uses SWIG to provide a python interface. 
The CRPropa functionality can be extended in Python in several ways.

### Implementing base classes
The base classes Module, Source and SourceProperty can be implemented in Python and used like the built-in classes. Here is an example for a custom Module:
```python
class MyModule(Module):
    """ Reduces the cosmic ray energy by 10% in each step """
    def process(self, c):
        c.current.setEnergy(c.current.getEnergy() * 0.9)

m = ModuleList()
m.add(MyModule())

c = Candidate()
c.current.setEnergy(10)
m.process(c)
print c.current.getEnergy()
```

When redefining the init function make sure to call the super classes init function as well.
```python
class MyModule(Module):
    def __init__(self):
        Module.__init__(self)
```


The initial properties of a cosmic rays can be set with a Source, composed of several SourceProperties.
Custom SourceProperties can be written in the following way:
```python
class MySourceProperty(SourceProperty):
    """ Set the initial energy to 10 EeV """
    def prepareParticle(self, particleState):
        particleState.setEnergy(10 * EeV)

s = Source()
s.addProperty(MySourceProperty())
c = s.getCandidate()
print c.current.getEnergy()
```

The redshift is stored in the Candidate, not in the ParticleState. To set it with a SourceProperty use the following:
```python
class MySourceProperty2(SourceProperty):
    """ Set the initial redshift """
    def prepareCandidate(self, candidate):
        candidate.setRedshift(0.6)

s = Source()
s.addProperty(MySourceProperty())
c = s.getCandidate()
print c.getRedshift()
```

### Replacing base classes
Alternatively the simulation chain (ModuleList) and the cosmic ray source (Source, SourceProperty) can be explicitly written down and replaced.
The following code replaces a ModuleList:
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
The created Candidates can either be propagated one-by-one, or first collected in a CandidateVector and then propagated.
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



 


