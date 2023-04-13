# Installation
## Download

Download and unzip the [latest release](https://github.com/CRPropa/CRPropa3/releases/latest) (recommended), or, alternatively, download the [current development snapshot](https://github.com/CRPropa/CRPropa3/archive/master.zip), or clone the repository with

```sh
git clone https://github.com/CRPropa/CRPropa3.git
```

## Prerequisites
+ C++ Compiler with C++11 support (gcc, clang and icc are known to work)
+ Fortran Compiler: to compile SOPHIA
+ numpy: for scientific computations

Optionally CRPropa can be compiled with the following dependencies to enable certain functionality.
+ Python and SWIG: to use CRPropa from python (tested for > Python 2.7 and > SWIG 3.0.4)
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
+ eigen: Linear algebra
+ healpix_base: Equal area pixelization of the sphere





## Build and Installation Variants
### Installation in system path

1. CRPropa uses CMAKE to configure the Makefile. From the build directory call
   ccmake or cmake. See the next section for a list of configuration flags.
    ```sh
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
    make
    make install
    ```

2. A set of unit tests can be run with ```make test```. If the tests are
   successful continue with ```make install``` to install CRPropa at the
   specified path, or leave it in the build directory.  Make sure the
   environment variables are set accordingly: E.g. for an installation under
   $HOME/.local and using Python 2.7 set
    ```sh
    export PATH=$HOME/.local/bin:$PATH
    export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages:$PYTHONPATH
    export PKG_CONFIG_PATH=$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH
    ```

However, we highly recommend to use a virtualenv setup to install CRPropa!


### Installation in python virtualenv
CRPropa is typically run on clusters where superuser access is not always
available to the user. Besides that, it is easier to ensure the reproducibility
of simulations in a user controlled and clean environment. Thus, the user
space deployment without privileged access to the system would be a preferred
way. Python provides the most flexible access to CRPropa features, hence,
Python and SWIG are required. To avoid clashes with the system's Python and its
libraries, Python virtual environment will be used as well.

This procedure brings a few extra steps compared to the already given plain
installation from source, but this kind of CRPropa deployment will be a
worthwhile effort afterwards.

1. Choose a location of the deployment and save it in an environment variable to avoid retyping, for example,
    ```sh
    export CRPROPA_DIR=$HOME"/.virtualenvs/crpropa"
    ```
    and make the directory
    ```sh
    mkdir -p $CRPROPA_DIR
    ```

2. Initialize the Python virtual environment with the virtualenv command,
    ```sh
    virtualenv $CRPROPA_DIR
    ```
    if there is virtualenv available on the system.
		If the virtualenv is not installed on a system, try to use your operating
		system software repository to install it (usually the package is called
		`virtualenv`, `python-virtualenv`, `python3-virtualenv` or
		`python2-virtualenv`). There is also an option to manually download it,
		un-zip it, and run it:
    ```sh
    wget https://github.com/pypa/virtualenv/archive/develop.zip
    unzip develop.zip
    python virtualenv-develop/virtualenv.py $CRPROPA_DIR
    ```

    Finally, activate the newly created virtual environment:
    ```sh
    source $CRPROPA_DIR"/bin/activate"
    ```

3. Check the dependencies and install at least mandatory ones (see [prerequisites](#prerequisites)). This can be done with package managers (see the [package list](#notes-for-specific-operating-systems) in different operating systems). If packages are installed from source, during the compilation the installation prefix should be specified:
    ```sh
    ./configure --prefix=$CRPROPA_DIR
    make
    make install
    ```

    To install python dependencies and libraries use `pip`. Example: `pip install numpy`.

4. Compile and install CRPropa (please note specific [insturctions for different operating systems](#notes-for-specific-operating-systems)).
    ```sh
    cd $CRPROPA_DIR
    git clone https://github.com/CRPropa/CRPropa3.git
    cd CRPropa3
    mkdir build
    cd build
    CMAKE_PREFIX_PATH=$CRPROPA_DIR cmake -DCMAKE_INSTALL_PREFIX=$CRPROPA_DIR ..
    make
    make install
    ```

5. A set of unit tests can be run with ```make test```. 

6. (optional) Check the installation.
    ```python
    python
    import crpropa
    ```
    The last command must execute without any output. To check if dependencies are installed and linked correctly use the following Python command, e.g. to test the availability of FFTW3:
    ```python
    'initTurbulence' in dir(crpropa)
    ```

There also exists [bash script](https://github.com/adundovi/CRPropa3-scripts/tree/master/deploy_crpropa) for GNU/Linux systems which automate the described procedure.


### CMake flags
When using cmake, the following options can be set by adding flags to the cmake command, e.g.
```
cmake -DENABLE_PYTHON=ON ..
```

+ Set the install path ```-DCMAKE_INSTALL_PREFIX=/my/install/path```
+ Enable Galactic magnetic lens ```-DENABLE_GALACTICMAGETICLENS=ON```
+ Enable FFTW3 (turbulent magnetic fields) ```-DENABLE_FFTW3F=ON```
+ Enable OpenMP (multi-core parallel computing) ```-DENABLE_OPENMP=ON```
+ Enable Python (Python interface with SWIG) ```-DENABLE_PYTHON=ON```
+ Enable HDF5 (HDF5 output) ```-DENABLE_HDF5=ON```
+ Enable [Quimby](https://git.rwth-aachen.de/3pia/forge/quimby) (multiresolution MHD fields) ```-DENABLE_QUIMBY=ON```
+ Enable the data file download (can be set to "off" if it is manually provided) ```-DDOWNLOAD_DATA=ON```
+ Enable unit-tests ```-DENABLE_TESTING=ON```
+ Enable Coverage (code coverage tool) ```-DENABLE_COVERAGE=ON```
+ Enable Git ```-DENABLE_GIT=ON```
+ Optimized parallelization usage for simulations with few particles ```-DOMP_SCHEDULE:STRING=dynamic``` (see [discussion](https://github.com/CRPropa/CRPropa3/issues/117))
+ Enable SWIG-builtin ```-DENABLE_SWIG_BUILTIN=ON```
+ Debugging symbols included: ```-DCMAKE_BUILD_TYPE:STRING=Debug```

  Generally, for compilers CMake recognise the following env variables: CC, CXX, FC. For example:
  ```
  export FC=/usr/bin/gfortran
  ```
  while CC and CXX are used C and C++ compilers, respectively.

+ Additional flags for Intel compiler
  ```
  -DCMAKE_SHARED_LINKER_FLAGS="-lifcore"
  -DCMAKE_Fortran_COMPILER=ifort
  ```

+ The PlaneWaveTurbulence computation can be improved using the FAST_WAVES flag (see [documentation](https://crpropa.github.io/CRPropa3/buildingblocks/MagneticFields.html#classcrpropa_1_1PlaneWaveTurbulence) for details):
  1. Enable FAST_WAVES flag ```-DFAST_WAVES=ON```
  2. Enable SIMD_EXTENSIONS ```-DSIMD_EXTENSIONS:STRING=native``` (the compiler will automatically detect support for your CPU and run the build with the appropriate settings).

  Note: If your CPU does not support the necessary extensions, the build will fail with an error telling you so. In this case, you won’t be able to use the optimization; go back into cmake, disable FAST_WAVES, and build again. If the build runs through without errors, the code is built with the optimization.

+ Quite often there are multiple Python versions installed in a system. This is likely the cause of many (if not most) of the installation problems related to Python. To prevent conflicts among them, one can explicitly refer to the Python version to be used. Example:
  ```
  -DCMAKE_PYTHON_EXECUTABLE=/usr/bin/python
  -DCMAKE_PYTHON_INCLUDE_DIR=<path_to_folder_containing_Python.h>
  -DCMAKE_PYTHON_LIBRARY=<path_to_file>/libpython<version_tag>.so
  ```
  Note that in systems running OSX, the extension .so should be replaced by .dylib.

## Notes for Specific Operating Systems

### Debian / Ubuntu
In a clean minimal **Ubuntu (17.10)** installation the following packages should be installed to build and run CRPropa with most of the options:
  ```sh
  sudo apt install python-virtualenv build-essential git cmake swig \
  gfortran python-dev fftw3-dev zlib1g-dev libmuparser-dev libhdf5-dev pkg-config
  ```

### Fedora/CentOS/RHEL
For Fedora/CentOS/RHEL the required packages to build CRPropa:
   ```sh
   yum install git cmake gcc gcc-gfortran gcc-c++ make swig zlib-devel \
   muParser-devel hdf5-devel fftw-devel python-devel
  ```
In case of CentOS/RHEL 7, the SWIG version is too old and has to be built from source.


### Mac OS X
Tested on version 12.5.1 with M1 pro where command line developer tools are installed. 
Install Python3, and llvm from Homebrew, and specify the following paths to the Python and llvm directories in the Homebrew folder after step 3 of the above installation, e.g. (please use your exact versions):
  ```sh
   export LLVM_DIR="/opt/homebrew/Cellar/llvm/15.0.7_1"
   PYTHON_VERSION=3.10
   LLVM_VERSION=15.0.7
   PYTHON_DIR=/opt/homebrew/Cellar/python@3.10/3.10.9/Frameworks/Python.framework/Versions/3.10
  ```
and replace the command in step 4 of the installation routine
  ```sh
  CMAKE_PREFIX_PATH=$CRPROPA_DIR cmake -DCMAKE_INSTALL_PREFIX=$CRPROPA_DIR ..
  ```
with
  ```sh
   cmake .. \
   -DCMAKE_INSTALL_PREFIX=$CRPROPA_DIR \
   -DPYTHON_EXECUTABLE=$PYTHON_DIR/bin/python$PYTHON_VERSION \
   -DPYTHON_LIBRARY=$PYTHON_DIR/lib/libpython$PYTHON_VERSION.dylib \
   -DPYTHON_INCLUDE_PATH=$PYTHON_DIR/include/python$PYTHON_VERSION \
   -DCMAKE_C_COMPILER=$LLVM_DIR/bin/clang \
   -DCMAKE_CXX_COMPILER=$LLVM_DIR/bin/clang++ \
   -DOpenMP_CXX_FLAGS="-fopenmp -I$LLVM_DIR/lib/clang/$LLVM_VERSION/include" \
   -DOpenMP_C_FLAGS="-fopenmp =libomp -I$LLVM_DIR/lib/clang/$LLVM_VERSION/include" \
   -DOpenMP_libomp_LIBRARY=$LLVM_DIR/lib/libomp.dylib \
   -DCMAKE_SHARED_LINKER_FLAGS="-L$LLVM_DIR/lib -lomp -Wl,-rpath,$LLVM_DIR/lib" \
   -DOpenMP_C_LIB_NAMES=libomp \
   -DOpenMP_CXX_LIB_NAMES=libomp \
   -DNO_TCMALLOC=TRUE
  ```
Check that all paths are set correctly with the following command in the build folder
  ```sh
   ccmake .. 
  ```
and configure and generate again after changes.


