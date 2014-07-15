CRPropa3
========

This is the development version of [CRPropa](https://crpropa.desy.de/Main_Page).
It features a very flexible setup of simulation, support for specialized extragalactic magnetic fields, galactic lensing, python - steering and parallelization.

## Install from source
1. Download the source
    ```
    git clone https://github.com/CRPropa/CRPropa3.git
    ```
2. Download the data repository
    ```
    git clone https://github.com/CRPropa/CRPropa3-data.git CRPropa3/data
    ```
3. CRPropa uses CMAKE to configure. From the build directory call ccmake or cmake. See the next section for a list of configuration flags
    ```
    cd build
    ccmake ..
    ```
4. After the configuration run make and make install as usual
    ```
    make
    make install
    ```

#### CMake flags
We recommend using ccmake to view and set the options through the user interface.
When using cmake, options can be set by adding flags to the cmake command, e.g. ```cmake -DENABLE_PYTHON=ON ..```

+ Set the install path
```-DCMAKE_INSTALL_PREFIX=/my/install/path```
+ Enable OpenMP
```-DENABLE_OPENMP=ON```
+ Enable Python
```-DENABLE_PYTHON=ON```
+ Enable ROOT
```-DENABLE_ROOT=ON```
+ Enable FFTW3
```-DENABLE_FFTW3F=ON```
+ Enable testing with googletest
```-DENABLE_TESTING=ON```
+ Additional flags for Intel compiler
```
-DCMAKE_SHARED_LINKER_FLAGS="-lifcore"
-DCMAKE_Fortran_COMPILER=ifort
```

#### Required dependencies
+ C++ Compiler
+ Fortran Compiler: to compile SOPHIA

#### Provided dependencies
+ SOPHIA: photo-hadronic interactions
+ EleCa and dint: electromagnetic cascades
+ googletest: unit-testing
+ HepPID: particle ID library
+ kiss: small tool collection
+ pugixml: for xml steering

#### Optional dependencies
+ Python and SWIG: to use CRPropa from python (tested for > Python 2.7 and > SWIG 2.0)
+ ROOT: for ROOT output
+ FFTW3: for turbulent magnetic field grids (FFTW3 with single precision is needed)
+ Gadget: magnetic fields for large scale structure data
+ OpenMP: for shared memory parallelization
+ googleperftools: for performance optimizations regarding shared memory parallelization
+ muparser: to define the source spectrum through a mathematical formula


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
    m.add(ElectronPairProduction(CMB))
    m.add(ElectronPairProduction(IRB))
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
