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

3. A set of unit tests can be run with ```make test```. If the tests are successful continue with ```make install``` to install CRPropa at the specified path, or leave it in the build directory.
Make sure the environment variables are set accordingly: E.g. for an installation under $HOME/.local and using Python 2.7 set
    ```sh
    export PATH=$HOME/.local/bin:$PATH
    export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages:$PYTHONPATH
    export PKG_CONFIG_PATH=$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH
    ```

#### Recommended set-up

There are different ways to install, set-up and use CRPropa, but for a common use case we would recommend the following one.

CRPropa is typically run on clusters where superuser access is not always available to the user. Besides that, it is easier to secure reproducibility of simulations in a user controlled and clean environment.
Thus, the deployment in a user space without privileged access to the system would be a preferred way. Python provides the most flexible access to CRPropa features, therefore, Python and SWIG will be required. To avoid clashes with the system's Python and its libraries, Python virtual environment will be used as well.

This procedure brings a few steps extra compared to the already given plain installation from the source, but later this kind of set-up will be a worthwhile effort.

1. Choose a location of the set-up and save it in the environment variable to avoid retyping:
`CRPROPA_DIR=$HOME"/.virtualenvs/crpropa"`
    and make the directory `mkdir -p $CRPROPA_DIR`

2. Initialize Python virtual environment with a virtualenv command.
    If the virtualenv is not already installed on a system (try `virtualenv` command), download it and un-zip it:
    ```sh
    wget https://github.com/pypa/virtualenv/archive/develop.zip
    unzip develop.zip
    python virtualenv-develop/virtualenv.py $CRPROPA_DIR
    ```

    Or instead of all this, just `virtualenv $CRPROPA_DIR` if there is virtualenv available on the system.
    
    Finally, activate the virtualenv:
    ```sh
    source $CRPROPA_DIR"/bin/activate"
    ```

3. Check for the dependencies and install required ones (see  [dependencies](#Dependencies) for details). During a compilation of the dependency given prefix is also needed:
    ```sh
    ./configure --prefix=$CRPROPA_DIR
    make
    make install
    ```
    To install python dependencies and libraries use `pip`. Example: `pip install numpy`.

4. Compile and install CRPropa.
    ```sh
    git clone https://github.com/CRPropa/CRPropa3.git
    cd CRPropa3
    mkdir build
    cd build
    CMAKE_PREFIX_PATH=$CRPROPA_DIR cmake -DCMAKE_INSTALL_PREFIX=$CRPROPA_DIR ..
    make
    make install
    ```

5. Add CRPropa to the virtualenv path.
    ```sh
    echo "export LD_LIBRARY_PATH=$CRPROPA_DIR/lib:\$LD_LIBRARY_PATH" >> $CRPROPA_DIR"/bin/activate"
    ```
    Then `deactivate` and activate virtualenv again: `source $CRPROPA_DIR"/bin/activate"`.

6. Check the set-up.
    ```python
    python
    import crpropa
    ```
    The last command must execute without any output. To check if dependencies are installed and linked correctly use the following Python command, e.g. to test the availability of FFTW3:
    ```python
    'initTurbulence' in dir(crpropa)
    ```

There also exists [bash script](https://github.com/adundovi/CRPropa3-scripts/tree/master/deploy_crpropa) for GNU/Linux systems which automate the described procedure.

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

#### <a name="Dependencies"></a>Dependencies
+ C++ Compiler
+ Fortran Compiler: to compile SOPHIA

Optionally CRPropa can be compiled with the following dependencies to enable certain functionality.
+ Python and SWIG: to use CRPropa from python (tested for > Python 2.7 and > SWIG 2.0.4)
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