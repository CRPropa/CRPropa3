# Installation

## Over Conda Package (recommended)

CRPropa can be installed simply over its [conda package](https://anaconda.org/channels/crpropa/packages/crpropa/overview) on every linux or macos system. For the most recent release simply use:

```sh
conda install -c conda-forge crpropa::crpropa
```

For the current master you can use:

```sh
conda install -c conda-forge crpropa::crpropa==master
```

You will find more information on how to install conda on their [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Testing CRPropa

You can test your CRPropa installation by just using:

```sh
testCRPropa
```

You can find the test logs at the following location:

```sh
$CONDA_PREFIX/share/crpropa/test/Testing/Temporary/
```

### Notes

The conda package does not come with all features:
- No local documentation
- No [coverage (lcov)](https://github.com/linux-test-project/lcov) report generation for tests
- No [QUIMBY](https://github.com/CRPropa/Quimby)
- No SIMD extensions for arm processors (not supported by the architecture)

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

## Building from Source

### Requirements
In general you need the following tools to build CRPropa3 from source, below are instructions how to install those on different operating systems.

- C++ Compiler with `C++17` support (`gcc` and `clang` are known to work)
- Fortran Compiler: to compile SOPHIA (for example `gfortran`)

Optionally CRPropa can be compiled with the following dependencies to enable certain functionality.
- `Python`, `NumPy`, and `SWIG`: to use CRPropa from python (tested for >= Python 3.10 and > SWIG 4.0.2)
- `FFTW3`: for turbulent magnetic field grids (FFTW3 with single precision is needed)
- `googleperftools`: for performance optimizations regarding shared memory parallelization
- `muparser`: to define the source spectrum through a mathematical formula
- `doxygen`: to build a `doxygen` documentation
- `lcov`, `genhtml`: to build coverage report with `cmake --build /path/to/your/buildfolder --target coverage` (requires executed tests over `ctest`)
- `sphinx`, `sphinx_rtd_theme`, `m2r2`, `nbsphinx`, `lxml_html_clean`, `breathe`, `pandoc`, `exhale`: to build this documentation from the `doxygen` generated documentation with `cmake --build /path/to/your/buildfolder --target doc` and possibly include the coverage report by copying the by `coverage` generated `coverageReport` to `doc/pages/coverageReport` and then do `cmake --build /path/to/your/buildfolder --target coverage`. You might want to install the mentioned packages over `pip` rather than `conda` since there is a known [bug](https://github.com/sphinx-doc/sphinx/issues/12239).
- `hdf5`: to enable the option to generate binary output

### CMake Flag Documentation

( name = default : explanation )
- `BUILD_DOC = OFF` : This enables the building of a `doxygen` version of the documentation, this will look very bare bone. To build a better documentation additionally install `sphinx` and use `cmake --build /path/to/your/buildfolder --target doc` while this options is set to `ON`.
- `CMAKE_INSTALL_PREFIX = /usr/local` : The installation prefix, a standard variable in every `cmake` build, this specifies where `cmake` should install the project. You should ensure you are actually allowed to write to that location.
- `DOWNLOAD_DATA = ON` : Whether to download the [data](https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz). You might want to disable this if you already have data downloaded or generated over [CRPrpa3-data](https://github.com/CRPropa/CRPropa3-data/tree/master) and do not want to override the existing data, or downloading and extracting is a problem for you.
- `ENABLE_COVERAGE = OFF` : Whether to enable [coverage (lcov)](https://github.com/linux-test-project/lcov) support. To build the coverage support you can use `cmake --build /path/to/your/buildfolder --target coverage` after doing `ctest`. You will also need `genhtml` and `lcov`
- `ENABLE_GIT = ON` : This just includes the corresponding CRPropa version correctly in the project
- `ENABLE_HDF5 = ON` : This enables the generation of HDF5 binary output, this does not replace the normal output.
- `ENABLE_OPENMP = ON` : Enables OpenMP parallelization, should cause a major speed up for any modern CPU.
- `ENABLE_PYTHON = ON` : Enables generation of python bindings over `SWIG`. If you want to use CRPropa fully from python you need to enable that option.
- `ENABLE_QUIMBY = OFF` : Enables [QUIMBY](https://github.com/CRPropa/Quimby) support, however you would need to build that tool yourself. You would want to include the quimby include path to your prefix path, in conda environments this is usually a given.
- `ENABLE_SWIG_BUILTIN = ON` : This enables us to create python-builtin types rather than proxies which increases performance.
- `ENABLE_TESTING = ON` : Enables the creation of test executables, the tests can then be done by using `ctest --output-on-failure`.
- `FAST_WAVES = OFF` : Enables the usage of SIMD extensions in `PlaneWaveTurbulence`, this can increase the performance by a lot if supported by your CPU. To check if you CPU supports this feature use `lscpu | grep -e avx -e fma`.
- `OMP_SCHEDULE = static,100` : The OMP strategy to use, to see more infos see [OMP Documentation](https://www.openmp.org/spec-html/5.0/openmpse49.html)
- `SIMD_EXTENSIONS = none` : The SIMD flag to use, allowed are `avx`, `avx+fma`, `native` and `none`. Check with `lscpu | grep -e avx -e fma` what is supported on your CPU, you could also use `native` to use whatever is available.
- `Python_INSTALL_PACKAGE_DIR` : With this variable you can specify your own install target for the python packages without changing the search path for the required sitelibs like `numpy`.

### Building on Linux

The following instructions are for a global installation, for a local installation see the conda section.
First you need to install the requirements with your favorite package manager, in the installation commands is every optional performance package included:

#### Debian/Ubuntu

```sh
sudo apt update
sudo apt install cmake cmake-curses-gui ninja-build git gcc g++ gfortran python3-dev python3-numpy-dev python-dev-is-python3 swig libhdf5-dev libmuparser-dev libfftw3-dev libgoogle-perftools-dev
```

#### Fedora

```sh
sudo dnf check-update
sudo dnf install cmake ccmake ninja git gcc g++ gfortran python3-devel python3-numpy swig hdf5-devel muParser-devel fftw3-devel gperftools
```

#### Arch

```sh
sudo pacman -Sy
sudo pacman -S cmake ninja git gcc gcc-fortran python python-numpy swig hdf5 muparser fftw gperftools pkgconfig
```

#### Conda

To see how to install conda please follow their [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

```sh
conda install -c conda-forge compilers git cmake ninja swig zlib gperftools fftw hdf5 muparser python numpy pkgconfig
```

You should also set `CMAKE_INSTALL_PREFIX=$CONDA_PREFIX` during the configure step.

#### Virtual Environment

If you want to avoid `conda` but still want to use additional python packages for example for the documentation,
you could also create a python virtual environment and install CRPropa over that.
To do that you need to do some additional steps, assuming you do not have set up a virtual environment yet:

```sh
sudo apt install python3-venv
python -m venv /path/to/your/venv
source /path/to/your/venv/bin/activate
pip install numpy
```

After these steps you can then just continue with the [building](#Building) step, the python paths pointing to the virtual environment should be prioritized.

In this example it is assumed you are on ubuntu for the installation of `python3-venv`, on fedora the corresponding package is called `python-venv` and on arch it is already included in the `python` package.

#### Building

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

# finally install; sudo might be needed
cmake --install .
```

### Building on MacOS

#### Brew

For a clean OS X (Sonoma 14+) installation, if you use [Homebrew](https://brew.sh/), the main dependencies can be installed as follows:

```sh
brew install cmake ninja gfortran libomp python swig numpy hdf5 fftw muparser gperftools
```

#### MacPorts

Similarly to Brew, if you use [MacPorts](https://www.macports.org/) instead of Homebrew, download the corresponding packages:

```sh
sudo port install cmake ninja gcc15 libomp python swig pkgconfig py-numpy hdf5 fftw-3 muparser gperftools
sudo port select --set gcc mp-gcc15
```

With MacPorts you might struggle to find the `python` and `numpy` libraries and headers, if so just follow [Notes](#notes-1).

#### Conda

In a conda environment you can build CRPropa3 [exactly like on a linux machine](#conda).

#### Building

With the now installed requirements we can now start building CRPropa. `libomp` is usually not found, so we need to specify the location and corresponding flags here:

```sh
# clone the repository wherever you want:
git clone https://github.com/CRPropa/CRPropa3.git

# build the project, you can change any of the above mentioned cmake flags with -D
# if you are not able to install ninja or problems arise use -G "Unix Makefiles"
cd CRPropa3
mkdir build && cd build
cmake .. -G Ninja \
  -DCMAKE_INSTALL_PREFIX="/path/to/your/installation" \
  -DOpenMP_CXX_FLAGS="-I/usr/local/opt/libomp/include" \
  -DOpenMP_C_FLAGS="-I/usr/local/opt/libomp/include" \
  -DOpenMP_libomp_LIBRARY=/usr/local/opt/libomp/lib/libomp.dylib \
  -DCMAKE_SHARED_LINKER_FLAGS="-L/usr/local/opt/libomp/lib -lomp -Wl,-rpath,/usr/local/opt/libomp/lib" \
  -DOpenMP_C_LIB_NAMES=libomp \
  -DOpenMP_CXX_LIB_NAMES=libomp
cmake --build . -j

# optionally do the tests to check if everything was build correctly
ctest --output-on-failure --repeat until-pass:3

# finally install; sudo might be needed
cmake --install .
```

### Building on Windows

Currently, we do not officially support Windows, it is advised to install UbuntuWSL and install it there.

## Notes

- Sometimes CMake has difficulties finding the correct `python` and `numpy` header and executables, to help CMake finding those define the following environment variables before using `cmake`:
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
  -DPython_INCLUDE_DIRS=${PYTHON_INCLUDE_DIR} \
  -DPython_NumPy_INCLUDE_DIRS=${NUMPY_INCLUDE_DIR} \
  -DPython_INSTALL_PACKAGE_DIR=${PYTHON_INSTALL_PACKAGE_DIR}
```
- To do the coverage report first install `lcov` and `genhtml`, then do the tests with `ctest` and built the coverage target with `cmake --build /path/to/your/buildfolder --target coverage`.
- It is generally advised to use `ninja` instead of the default `make`, use it with the `-G Ninja` flag for `cmake`, you can install it over `pip`, `conda` or your package manager
- You might want to generate some python stubs to get code recommendation from tools like `pylance`, for that install `pybind11-stubgen` over `pip`, `conda` or your package manager and then do:
```sh
pybind11-stubgen -o $(python -c "import sysconfig; print(sysconfig.get_paths()['purelib'])") crpropa
```

# Known Issues

- The command `testCRPropa` is based on `ctest` which is provided by `cmake` which is a runtime requirement for CRPropa for this reason. However, since we do not require a specific version, it is possible, that the provided `ctest` is too modern for your machine and runs into issues even though `crpropa` might still be runnable.
The current workaround would be either to install an older `ctest` version or to run the tests manually which are located in `$CONDA_PREFIX/share/crpropa/test`.