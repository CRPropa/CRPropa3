CRPropa3
========

This is the development version of [CRPropa](https://crpropa.desy.de/Main_Page).
It features a very flexible setup of simulation, support for specialized extragalactic magnetic fields, galactic lensing, python - steering and parallelization.

## Install from source
1. Download the data repository
```
cd CRPropa3
git clone https://github.com/CRPropa/CRPropa3-data.git data
```
2. CRPropa uses CMAKE to configure. From the build directory call ccmake or cmake
```
cd build
ccmake ..
```
3. Afterward configuring run make and make install as usual
```
make
make install
```

The install path can be set with -DCMAKE_INSTALL_PREFIX=/my/path or with the option browser when using ccmake.

Notes for Intel Compiler:
use -DCMAKE_SHARED_LINKER_FLAGS="-lifcore" -DCMAKE_Fortran_COMPILER=ifort

#### Provided dependencies
+ SOPHIA
    + for photo-hadronic interactions
+ googletest
    + for unit-tests

#### Optional dependencies
+ Python and SWIG
    + to use CRPropa from Python
    + tested for > Python 2.7
    + tested for > SWIG 2.0
+ FFTW3F
    + for turbulent magnetic field grids
    + CRPropa needs the FFTW3 library compiled with the single precision option
+ Gadget
    + Magnetic fields for large scale structure data
+ OpenMP
    + for shared memory parallelization
+ googleperftools
    + for performance optimizations regarding shared memory parallelization


## Getting started
We recommend using CRPropa via python. The documentation can be found on the [CRPropa wiki](https://crpropa.desy.de/CRPropa3]).
For a 1D simulation try

 ```python
    from crpropa import *

    # simulation setup
    m = ModuleList()
    m.add(SimplePropagation(0, 10 * Mpc))
    m.add(PhotoPionProduction(CMB))
    m.add(PhotoPionProduction(IRB))
    m.add(PhotoDisintegration(CMB))
    m.add(PhotoDisintegration(IRB))
    m.add(ElectronPairProduction(CMB_IRB))
    m.add(NuclearDecay())
    m.add(MinimumEnergy(1 * EeV))
    m.add(Observer1D())
    m.add(EventOutput1D('events.txt'))

    # source setup
    source = Source()
    source.addProperty(SourceParticleType(nucleusId(56, 26)))
    source.addProperty(SourcePowerLawSpectrum(1 * EeV, 500 * EeV, -2))
    source.addProperty(SourceUniform1D(3 * Mpc, 2000 * Mpc))

    # run simulation
    m.setShowProgress(True)
    m.run(source, 1000, True)
 ```

## CRPropa 2 compatibility
For backwards compatibility CRPropa 3 can be steered via XML cards ($cropra-xmlrun some_steeringcard.xml)
However, CRPropa 2 does not fully enforce the XML-1.0 standard (http://www.w3.org/XML).
To comply with the standard a few modifications to exisisting steering cards might have to be made.
Modification are

1. XML-cards can have only one root node. The root node is

        <CRPropa>
        ...
        </CRPropa>

2. All nodes including the optional header node need to be closed, e.g.

        <?xml version="1.0" standalone=no?>
        <Option1> ... </Option1> or
        <Option2/>

3. Values need to be embraced by quotes, e.g.

        <NumTrajectories="1000"/>
