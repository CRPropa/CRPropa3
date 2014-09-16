#### Install from source
1. Download and unzip the [source file](https://github.com/CRPropa/CRPropa3/archive/master.zip) or clone the repository with
    ```
    git clone https://github.com/CRPropa/CRPropa3.git
    ```

2. CRPropa uses CMAKE to configure the Makefile. From the build directory call ccmake or cmake. See the next section for a list of configuration flags.
    ```
    mkdir build
    cd build
    ccmake ..
    make
    ```

  CMake also automatically downloads a tarball with the interaction data. If the download fails, please manually download the specified file, rename it to data.tar.gz and place it in the source folder before continuing with cmake.

3. A set of unit tests can be run with ```make test```. If the tests are successful continue with ```make install``` to install CRPropa at the specified path, or leave it in the build directory.
Make sure the environment variables are set accordingly: E.g. for an installation under $HOME/.local and using Python 2.7 set
    ```sh
    export PATH=$HOME/.local/bin:$PATH
    export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages:$PYTHONPATH
    export PKG_CONFIG_PATH=$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH
    ```

#### CMake flags
We recommend using ccmake to view and set the options through the user interface.
When using cmake, options can be set by adding flags to the cmake command, e.g. 
```
cmake -DENABLE_PYTHON=ON ..
```

+ Set the install path ```-DCMAKE_INSTALL_PREFIX=/my/install/path```
+ Enable OpenMP ```-DENABLE_OPENMP=ON```
+ Enable Python ```-DENABLE_PYTHON=ON```
+ Enable ROOT ```-DENABLE_ROOT=ON```
+ Enable FFTW3 ```-DENABLE_FFTW3F=ON```
+ Enable unit-tests ```-DENABLE_TESTING=ON```
+ Additional flags for Intel compiler
  ```
  -DCMAKE_SHARED_LINKER_FLAGS="-lifcore"
  -DCMAKE_Fortran_COMPILER=ifort
  ```

#### Dependencies
+ C++ Compiler
+ Fortran Compiler: to compile SOPHIA

Optionally CRPropa can be compiled with the following dependencies to enable certain functionality.
+ Python and SWIG: to use CRPropa from python (tested for > Python 2.7 and > SWIG 2.0)
+ ROOT: for ROOT output
+ FFTW3: for turbulent magnetic field grids (FFTW3 with single precision is needed)
+ Gadget: magnetic fields for large scale structure data
+ OpenMP: for shared memory parallelization
+ googleperftools: for performance optimizations regarding shared memory parallelization
+ muparser: to define the source spectrum through a mathematical formula

The following packages are provided with the source code and do not need to be installed separately.
+ SOPHIA: photo-hadronic interactions
+ EleCa and dint: electromagnetic cascades
+ googletest: unit-testing
+ HepPID: particle ID library
+ kiss: small tool collection
+ pugixml: for xml steering