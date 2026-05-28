# Installation

## Over Conda Package

CRPropa can be installed simply over its [conda package]() on every linux system. For the most recent release simply use:

```sh
conda install crpropa::crpropa
```

For the current master you can use:

```sh
conda install crpropa::crpropa==master
```

### Testing CRPropa

You can test your CRPropa installation by just using:

```sh
testCRPropa
```

### Notes

The conda package does not come with all possible available features:
- No local documentation
- No [coverage (lcov)](https://github.com/linux-test-project/lcov) report generation for tests
- No [QUIMBY](https://github.com/CRPropa/Quimby)
- No SIMD extensions for arm processors (not supported by the achitecture)

If you want to include your custom data the location of the different data folders and files is

```sh
$CONDA_PREFIX/share/crpropa
```

If you need the location of the swig headers, they can be found at:

```sh
$CONDA_PREFIX/share/crpropa/swig_interface
```

the libraries and headers can be found at

```sh
$CONDA_PREFIX/lib
$CONDA_PREFIX/include
```

# Building from Source

## Requirements
In general you need the following tools to build CRPropa3 from source, below are instructions how to install those on different operating systems.

- C++ Compiler with `C++17` support (`gcc`, `clang` and `icc` are known to work)
- Fortran Compiler: to compile SOPHIA (for example `gfortran`)

Optionally CRPropa can be compiled with the following dependencies to enable certain functionality.
- `Python`, `NumPy`, and `SWIG`: to use CRPropa from python (tested for >= Python 3.10 and > SWIG 4.0.2)
- `FFTW3`: for turbulent magnetic field grids (FFTW3 with single precision is needed)
- `googleperftools`: for performance optimizations regarding shared memory parallelization
- `muparser`: to define the source spectrum through a mathematical formula
- `doxygen`: to build a `doxygen` documentation
- `lcov`, `genhtml`: to build coverage report with `cmake --build /path/to/your/buildfolder --target coverage` (requires executed tests over `ctest`)
- `sphinx`: to build this documentation from the `doxygen` generated documentation and possibly include the coverage report by copying the by `coverage` generated `coverageReport` to `doc/pages/coverageReport` and then do `cmake --build /path/to/your/buildfolder --target coverage`
- `hdf5`: to enable the option to generate binary output

## CMake Flag Documentation

( name = default : explanation )
- `BUILD_DOC = OFF` : This enables the building of a `doxygen` version of the documentation, this will look very bare bone. To build a better documentation additionally install `sphinx` and use `cmake --build /path/to/your/buildfolder --target doc` while this options is set to `ON`.
- `CMAKE_INSTALL_PREFIX = /usr/local` : The installation prefix, a standard variable in every `cmake` build, this specifies where `cmake` should install the project. You should ensure you are actually allowed to write to that location.
- `DOWNLOAD_DATA = ON` : Whether to download the [data](https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz). You might want to disable this if you already have data downloaded or generated over [CRPrpa3-data](https://github.com/CRPropa/CRPropa3-data/tree/master) and do not want to override the existing data, or downloading and extracting is a problem for you.
- `EIGEN_PATH = ""` : Here you can specify the path to your own eigen headers. However, CRPropa already provides a `eigen` version.
- `ENABLE_COVERAGE = OFF` : Whether to enable [coverage (lcov)](https://github.com/linux-test-project/lcov) support. To build the coverage support you can use `cmake --build /path/to/your/buildfolder --target coverage` after doing `ctest`. You will also need `genhtml` and `lcov`
- `ENABLE_GIT = ON` : This just includes the corresponding CRPropa version correctly in the project
- `ENABLE_HDF5 = ON` : This enables the generation of HDF5 binary output, this does not replace the normal output.
- `ENABLE_OPENMP = ON` : Enables OpenMP parallelization, should cause a major speed up for any modern CPU.
- `ENABLE_PYTHON = ON` : Enables generation of python bindings over `SWIG`. If you want to use CRPropa fully from python you need to enable that option.
- `ENABLE_QUIMBY = OFF` : Enables [QUIMBY](https://github.com/CRPropa/Quimby) support, however you would need to build that tool yourself. You would want to include the quimby include path to your prefix path, in conda environments this is usually a given.
- `ENABLE_SWIG_BUILTIN = ON` : This enables us to create python-builtin types rather than proxies which increases performance.
- `ENABLE_TESTING = ON` : Enables the creation of test executables, the tests can then be done by using `ctest --output-on-failure`.
- `FAST_WAVES = OFF` : Enables the usage of SIMD extensions in `PaneWaveTurbulence`, this can increase the performance by a lot if supported by your CPU. To check if you CPU supports this feature use `lscpu | grep -e avx -e fma`.
- `INSTALL_EIGEN = OFF` : Whether to install the provided eigen or not, this might override any preexisting version.
- `OMP_SCHEDULE = static,100` : The OMP strategy to use, to see more infos see [OMP Documentation](https://www.openmp.org/spec-html/5.0/openmpse49.html)
- `SIMD_EXTENSIONS = none` : The SIMD flag to use, allowed are `avx`, `avx+fma`, `native` and `none`. Check with `lscpu | grep -e avx -e fma` what is supported on your CPU, you could also use `native` to use whatever is available.

## Building on Linux

The following instructions are for a global installation, for a local installation see the conda section.
First you need to install the requirements with your favorite package manager, in the installation commands is every optional performance package included:

### Debian/Ubuntu

```sh
sudo apt update
sudo apt install cmake cmake-curses-gui ninja-build git gcc g++ gfortran python3-dev python3-numpy-dev python-dev-is-python3 swig libhdf5-dev libmuparser-dev libfftw3-dev libgoogle-perftools-dev
```

### Fedora

```sh
sudo dnf makecache
sudo dnf install cmake ccmake ninja git gcc g++ gfortran python3-devel python3-numpy swig hdf5-devel muParser-devel fftw3-devel gperftools
```

### Arch

```sh
sudo pacman -Sy
sudo pacman -S cmake ninja git gcc gcc-fortran python python-numpy swig hdf5 muparser fftw gperftools pkgconfig
```

### Conda

To see how to install conda please follow their [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

```sh
conda install -c conda-forge compilers git cmake ninja swig zlib gperftools fftw hdf5 muparser eigen python numpy pkgconfig
```

You should also set `CMAKE_INSTALL_PREFIX=$CONDA_PREFIX` during the configure step.

### Building

To then actually build the project do (independent of distribution):

```sh
# clone the repository wherever you want:
git clone https://github.com/CRPropa/CRPropa3.git

# build the project, you can change any of the above mentioned cmake flags with -D
# if you are not able to install ninja or problems arise use -G "Unix Makefiles"
cd CRPropa3
mkdir build && cd build
cmake .. -G Ninja -DCMAKE_INSTALL_PREFIX=/path/to/your/desired/install/location
cmake --build . -j

# optionally do the tests to check if everything was build correctly
ctest --output-on-failure --repeat until-pass:3

# finally install, if you did not specify a CMAKE_INSTALL_PREFIX where you have write rights
# of are not in a venv or conda environment you need to use the next line with sudo:
cmake --install .
```

## Building on MacOS

The following constructions are taken from an old version of this installation guide, the instructions could not be tested due to missing hardware.

For a clean OS X (Sonoma 14+) installation, if you use Homebrew, the main dependencies can be installed as follows:
   ```sh
   brew install hdf5 fftw cfitsio muparser libomp numpy swig
  ```
Similarly, if you use MacPorts instead of Homebrew, download the corresponding packages:
   ```sh
   sudo port install hdf5 fftw cfitsio muparser libomp numpy swig
  ```
Note that if you are using a Mac with Arm64 architecture (M1, M2, or M3 processors), `SIMD_EXTENSIONS` might not run straight away.


Some combinations of versions of the Apple's clang compiler and python might lead to installation errors.
In these cases, the user might want to consider the workaround below (tested on version 12.5.1 with M1 pro where command line developer tools are installed).

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
   -DPython_EXECUTABLE=$PYTHON_DIR/bin/python$PYTHON_VERSION \
   -DPython_LIBRARY=$PYTHON_DIR/lib/libpython$PYTHON_VERSION.dylib \
   -DPython_INCLUDE_PATH=$PYTHON_DIR/include/python$PYTHON_VERSION \
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

## Building on Windows

Currently, we do not officially support Windows, it is advised to install UbuntuWSL and install it there.
However, it might be possible to install it natively on Windows, you would however need to change some includes
that are only available on linux to the corresponding packages on Windows. It probably helps to use conda as a package manager.

## Notes

- Sometimes CMake has difficulties finding the correct python and numpy header and executables, to help CMake finding those define the following environment variables before using `cmake`:
```sh
export PYTHON_EXECUTABLE=$(which python)
export PYTHON_INCLUDE_DIR=$(${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_paths()['include'])")
export NUMPY_INCLUDE_DIR=$(${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())")
export PYTHON_INSTALL_PACKAGE_DIR=$(${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_paths()['purelib'])")
```
And hand them over to `cmake` as CMake-Variables so `find_python` can recognize them:
```sh
cmake .. \
	-DPython_EXECUTABLE=${PYTHON_EXECUTABLE} \
	-DPython_INCLUDE_DIR=${PYTHON_INCLUDE_DIR} \
	-DPython_NumPy_INCLUDE_DIR=${NUMPY_INCLUDE_DIR} \
	-DPython_INSTALL_PACKAGE_DIR=${PYTHON_INSTALL_PACKAGE_DIR}
```
- To do the coverage report first install `lcov` and `genhtml`, then do the tests with `ctest` and built the coverage target with `cmake --build /path/to/your/buildfolder --target coverage`.
- It is generally advised to use `ninja` instead of the default `make`, use it with the `-G Ninja` flag for `cmake`, you can install it over `pip`, `conda` or your package manager
- You might want to generate some python stubs to get code recommendation from tools like `pylance`, for that install `pybind11-stubgen` over `pip`, `conda` or your package manager and then do:
```sh
pybind11-stubgen -o $(python -c "import sysconfig; print(sysconfig.get_paths()['purelib'])") crpropa
```