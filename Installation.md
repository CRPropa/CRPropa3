### <a name="Download">Download the source code

Download and unzip the [latest release](https://github.com/CRPropa/CRPropa3/releases/latest) or the
[current development snapshot](https://github.com/CRPropa/CRPropa3/archive/master.zip), or clone the
repository with

```sh
git clone https://github.com/CRPropa/CRPropa3.git
```

### <a name="Install">Install from source (experienced users)

1. CRPropa uses CMAKE to configure the Makefile. From the build directory call
   ccmake or cmake. See the next section for a list of configuration flags.
    ```sh
    mkdir build
    cd build
    ccmake ..
    make
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

We highly recommend to use a virtualenv setup to install CRPropa!

### <a name="virtualenv">Install from source in python virtualenv - a detailed guide (RECOMMENDED)

CRPropa is typically run on clusters where superuser access is not always available to the user. Besides that, it is easier to ensure the reproducibility of simulations in a user controlled and clean environment.
Thus, the user space deployment without privileged access to the system would be a preferred way. Python provides the most flexible access to CRPropa features, hence, Python and SWIG are required. To avoid clashes with the system's Python and its libraries, Python virtual environment will be used as well.

This procedure brings a few extra steps compared to the already given plain installation from source, but this kind of CRPropa deployment will be a worthwhile effort afterwards.

1. Choose a location of the deployment and save it in an environment variable to avoid retyping, for example, 
    ```sh
    export CRPROPA_DIR=$HOME"/.virtualenvs/crpropa"
    ```
    and make the directory
    ```sh
    mkdir -p $CRPROPA_DIR`
    ```

2. Initialize the Python virtual environment with the virtualenv command.
    ```sh
    virtualenv $CRPROPA_DIR`
    ```
    if there is virtualenv available on the system.
    If the virtualenv is not installed on a system, try to use your operating system software repository to install it (usually the package is called `virtualenv`, `python-virtualenv`, `python3-virtualenv` or `python2-virtualenv`). There is also an option to manually download it, un-zip it and run it:
    ```sh
    wget https://github.com/pypa/virtualenv/archive/develop.zip
    unzip develop.zip
    python virtualenv-develop/virtualenv.py $CRPROPA_DIR
    ```

    Finally, activate the newly created virtual environment:
    ```sh
    source $CRPROPA_DIR"/bin/activate"
    ```

3. Check the dependencies (see  [dependencies](#Dependencies) for details) and install at least mandatory ones. This can be done with package managers (see the [package list](#dependencies-in-different-oses) in different OSes). If packages are installed from source, during the compilation the installation prefix should be specified:
    ```sh
    ./configure --prefix=$CRPROPA_DIR
    make
    make install
    ```
    
    To install python dependencies and libraries use `pip`. Example: `pip install numpy`.

4. Compile and install CRPropa.
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

5. (optional) Check the installation.
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

### <a name="Dependencies"></a>Dependencies
+ C++ Compiler (gcc, clang and icc are known to work)
+ Fortran Compiler: to compile SOPHIA

Optionally CRPropa can be compiled with the following dependencies to enable certain functionality.
+ Python and SWIG: to use CRPropa from python (tested for > Python 2.7 and > SWIG 2.0.4)
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

#### <a name="OS"></a>Dependencies in different OSes:

In a clean minimal **Ubuntu (17.10)** installation the following packages should be installed to build and run CRPropa with most of the options:
  ```sh
  sudo apt install python-virtualenv build-essential git cmake swig \
  gfortran python-dev fftw3-dev zlib1g-dev libmuparser-dev libhdf5-dev pkg-config
  ```

