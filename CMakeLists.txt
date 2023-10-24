if(APPLE)
  # rpath specific patches
  cmake_minimum_required(VERSION 2.8.12)
else(APPLE)
  # require > 2.8.8 for FILE DOWNLOAD fixes
  # allow < 2.8.12 for debian backports
  cmake_minimum_required(VERSION 2.8.11)
endif(APPLE)

project(CRPropa Fortran C CXX)
set(CRPROPA_RELEASE_VERSION 3.2+)  # Update for new release

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CRPROPA_EXTRA_SOURCES)
set(CRPROPA_EXTRA_INCLUDES)
set(CRPROPA_EXTRA_LIBRARIES)
set(CRPROPA_SWIG_DEFINES)
set(CRPROPA_SWIG_INPUTS)

macro(USE_CXX11)
  if(CMAKE_VERSION VERSION_LESS "3.1")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
    endif()
  else()
    set(CMAKE_CXX_STANDARD 11)
  endif()
endmacro(USE_CXX11)
USE_CXX11()

if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
  message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)

SET(SIMD_EXTENSIONS "none" CACHE STRING "Choose which of the SIMD instruction set extensions your target CPU supports. Possible values are \"native\" (use everything that's supported by the CPU you're building on), \"none\", \"avx\", and \"avx+fma\".")

if(SIMD_EXTENSIONS STREQUAL "none")
  message("With SIMD_EXTENSIONS \"none\" tricubic interpolation of vector grids is not possible. You should set SIMD_EXTENSION to a compatible value (\"avx\", \"avx+fma\", or -- depending on the build CPU -- \"native\"). Trilinear interpolation is still possible.")
else()
  add_definitions(-DHAVE_SIMD)
endif()

if(SIMD_EXTENSIONS STREQUAL "avx")
  SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -msse4.1 -msse4.2 -mavx" )
elseif(SIMD_EXTENSIONS STREQUAL "avx+fma")
  SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -msse4.1 -msse4.2 -mavx -mfma" )
elseif(SIMD_EXTENSIONS STREQUAL "native")
  SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -march=native" )
elseif(NOT SIMD_EXTENSIONS STREQUAL "none")
  message(SEND_ERROR "SIMD_EXTENSIONS must have one of these values: \"native\", \"none\", \"avx\", or \"avx+fma\".")
endif()

SET(FAST_WAVES OFF CACHE BOOL "Enable SIMD optimizations for PlaneWaveTurbulence. Requires that SIMD_EXTENSIONS is set to a compatible value (\"avx\", \"avx+fma\", or -- depending on the build CPU -- \"native\"). The build will fail if this option is set but SIMD_EXTENSIONS is not compatible, so you can adjust your settings.")
if(FAST_WAVES)
  if(SIMD_EXTENSIONS STREQUAL "none")
    message(SEND_ERROR "Enabling FAST_WAVES while SIMD_EXTENSIONS is set to \"none\" won't work, since SIMD is required for FAST_WAVES. Please try disabling FAST_WAVES, or try using another option for SIMD_EXTENSIONS.")
  else()
    add_definitions(-DFAST_WAVES)
  endif()
endif(FAST_WAVES)

# Add build type for profiling
SET(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -ggdb -fno-omit-frame-pointer")
# Enable extra warnings on debug builds
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wextra")

# Set default build-type to release to enable performance improvements
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")


# ----------------------------------------------------------------------------
# Version info from Git
# ----------------------------------------------------------------------------
# Embed the information about CRPropa version into the binaries for easier
# debugging with a reproducibility in mind
option(ENABLE_GIT "Embedding information about CRPropa version from git" ON)

if(ENABLE_GIT)
  include(GetGitRevisionDescription)
  get_git_head_revision(GIT_REFSPEC GIT_SHA1)
  git_describe(GIT_DESC --tags)
endif(ENABLE_GIT)
IF(NOT ENABLE_GIT OR (${GIT_SHA1} STREQUAL "GITDIR-NOTFOUND"))
  set(GIT_REFSPEC "")
  set(GIT_SHA1 "")
  set(GIT_DESC "${CRPROPA_RELEASE_VERSION}-no-git")
endif()

message(STATUS "CRPropa version: ${GIT_DESC} ${GIT_SHA1} ${GIT_REFSPEC}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/Version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/src/Version.cpp" @ONLY)
list(APPEND CRPROPA_EXTRA_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/src/Version.cpp")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# googletest (provided, see https://code.google.com/p/googletest/wiki/FAQ
#     Why is it not recommended use a pre-compiled copy of Google Test?)
option(ENABLE_TESTING "Build tests and enable test target" ON)
if(ENABLE_TESTING)
  include_directories(libs/gtest/include)
  add_subdirectory(libs/gtest)
  if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
  endif(APPLE)
endif(ENABLE_TESTING)

#
OPTION(ENABLE_COVERAGE "Add code for coverage analysis." OFF)
if(ENABLE_COVERAGE)
  find_program(LCOV_PATH lcov)
  find_program(GENHTML_PATH genhtml)
  if(LCOV_PATH AND GENHTML_PATH)
    message("Enabling coverage report via $make test && make coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    list(APPEND CRPROPA_EXTRA_LIBRARIES "-lgcov")
    list(APPEND CRPROPA_EXTRA_LIBRARIES "-fprofile-arcs")

    if(ENABLE_TESTING)
      SET(COVERAGE_LIBS "-coverage -lgcov")
      add_custom_target(coverage_clean
        COMMAND ${LCOV_PATH} --directory . --zerocounters
      )
      add_custom_target(coverage
        # generate coverage data
        COMMAND ${LCOV_PATH} --directory . --capture --output-file coverage.info VERBATIM
        # clean external libs
        COMMAND ${LCOV_PATH} --remove coverage.info "/usr/include/*" "/usr/lib/*" "*/libs/gtest/*" "*/libs/eigen3/*" "*/libs/zstream-cpp/*" "*/build/*" -o coverage.info.cleaned VERBATIM
        # Generate html output
        COMMAND ${GENHTML_PATH} -o coverageReport coverage.info.cleaned VERBATIM
        COMMAND echo "Generated coverage report in coverageReport/index.html"
      )
    endif(ENABLE_TESTING)

  else(LCOV_PATH AND GENHTML_PATH)
    if(NOT LCOV_PATH)
      message(WARNING "lcov not found, coverage report generation  not possible!")
    endif(NOT LCOV_PATH)
    if(NOT GENHTML_PATH)
      message(WARNING "genhtml not found, coverage report generation not possible!")
    endif(NOT GENHTML_PATH)
  endif(LCOV_PATH AND GENHTML_PATH)
endif(ENABLE_COVERAGE)

# kiss (provided)
add_subdirectory(libs/kiss)
list(APPEND CRPROPA_EXTRA_LIBRARIES kiss)
list(APPEND CRPROPA_EXTRA_INCLUDES libs/kiss/include)

# HepID (provided)
add_subdirectory(libs/HepPID)
list(APPEND CRPROPA_EXTRA_LIBRARIES HepPID)
list(APPEND CRPROPA_EXTRA_INCLUDES libs/HepPID/include)

# SOPHIA (provided)
add_subdirectory(libs/sophia)
list(APPEND CRPROPA_EXTRA_LIBRARIES sophia gfortran)
list(APPEND CRPROPA_EXTRA_INCLUDES libs/sophia)

# DINT (provided)
add_subdirectory(libs/dint)
list(APPEND CRPROPA_EXTRA_LIBRARIES dint)
list(APPEND CRPROPA_EXTRA_INCLUDES libs/dint/include)

# EleCa (provided)
add_subdirectory(libs/EleCa)
list(APPEND CRPROPA_EXTRA_LIBRARIES eleca)
list(APPEND CRPROPA_EXTRA_INCLUDES libs/EleCa/include)

# GlacticMagneticLenses
option(ENABLE_GALACTICMAGETICLENS "Galactic Magnetic Lens" ON)
option(INSTALL_EIGEN "Install provided EIGEN headers" OFF)
SET(EIGEN_PATH  "" CACHE STRING "Use EIGEN from this path instead of the version shipped with CRPropa")
SET(WITH_GALACTIC_LENSES FALSE)
if(ENABLE_GALACTICMAGETICLENS)
  SET(WITH_GALACTIC_LENSES TRUE)

  if(EIGEN_PATH)
    message("Using EIGEN from " ${EIGEN_PATH})
    list(APPEND CRPROPA_EXTRA_INCLUDES ${EIGEN_PATH})
  else(EIGEN_PATH)
    # Eigen redux (provided)
    message("Using provided EIGEN")
    list(APPEND CRPROPA_EXTRA_INCLUDES libs/eigen3)
  endif(EIGEN_PATH)
  if(INSTALL_EIGEN)
    install(DIRECTORY libs/eigen3/ DESTINATION include)
  endif(INSTALL_EIGEN)

  # healpix redux (provided)
  add_subdirectory(libs/healpix_base)
  list(APPEND CRPROPA_EXTRA_LIBRARIES healpix_base)
  list(APPEND CRPROPA_EXTRA_INCLUDES libs/healpix_base/include)
  install(DIRECTORY libs/healpix_base/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")

  list(APPEND CRPROPA_SWIG_DEFINES -DWITH_GALACTIC_LENSES)

  list(APPEND CRPROPA_EXTRA_SOURCES src/magneticLens/MagneticLens.cpp)
  list(APPEND CRPROPA_EXTRA_SOURCES src/magneticLens/ModelMatrix.cpp)
  list(APPEND CRPROPA_EXTRA_SOURCES src/magneticLens/Pixelization.cpp)
  list(APPEND CRPROPA_EXTRA_SOURCES src/magneticLens/ParticleMapsContainer.cpp)
endif(ENABLE_GALACTICMAGETICLENS)

# OpenMP (optional for shared memory multiprocessing)
option(ENABLE_OPENMP "OpenMP for multithreading" ON)
if(ENABLE_OPENMP)
  include(FindOpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif(OPENMP_FOUND)
endif(ENABLE_OPENMP)

# Additional configuration OMP_SCHEDULE
set(OMP_SCHEDULE "static,100" CACHE STRING "FORMAT type,chunksize")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/ModuleList.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/src/ModuleList.cpp" @ONLY)
list(APPEND CRPROPA_EXTRA_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/src/ModuleList.cpp")

# Google Performance Tools (optional as possible performance tweak for OpenMP)
find_package(GooglePerfTools)
set(TCMALLOC)
if(GOOGLE_PERFTOOLS_FOUND AND NOT NO_TCMALLOC)
  set(TCMALLOC ${TCMALLOC_LIBRARY})
  list(APPEND CRPROPA_EXTRA_LIBRARIES profiler)
endif(GOOGLE_PERFTOOLS_FOUND AND NOT NO_TCMALLOC)

# FFTW3F (optional for turbulent magnetic fields)
find_package(FFTW3F)
if(FFTW3F_FOUND)
  list(APPEND CRPROPA_EXTRA_INCLUDES ${FFTW3F_INCLUDE_DIR})
  list(APPEND CRPROPA_EXTRA_LIBRARIES ${FFTW3F_LIBRARY})
  add_definitions(-DCRPROPA_HAVE_FFTW3F)
  list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_FFTW3F)
endif(FFTW3F_FOUND)

# Quimby (optional for SPH magnetic fields)
option(ENABLE_QUIMBY "Quimby Support" ON)
if(ENABLE_QUIMBY)
  find_package(Quimby)
  if(QUIMBY_FOUND)
    list(APPEND CRPROPA_EXTRA_INCLUDES ${QUIMBY_INCLUDE_DIR})
    list(APPEND CRPROPA_EXTRA_LIBRARIES ${QUIMBY_LIBRARY})
    add_definitions (-DCRPROPA_HAVE_QUIMBY)
    list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_QUIMBY)
    list(APPEND SWIG_INCLUDE_DIRECTORIES ${QUIMBY_INCLUDE_DIR}/../share/quimby)
    list(APPEND SWIG_INCLUDE_DIRECTORIES ${QUIMBY_INCLUDE_DIR})
  endif(QUIMBY_FOUND)
endif(ENABLE_QUIMBY)

# SAGA (optinal for AMR grids)
find_package(SAGA)
if(SAGA_FOUND)
  list(APPEND CRPROPA_EXTRA_INCLUDES ${SAGA_INCLUDE_DIR})
  list(APPEND CRPROPA_EXTRA_LIBRARIES ${SAGA_LIBRARY})
  add_definitions (-DCRPROPA_HAVE_SAGA)
  list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_SAGA)
  list(APPEND SWIG_INCLUDE_DIRECTORIES ${SAGA_INCLUDE_DIR}/../python)
  list(APPEND SWIG_INCLUDE_DIRECTORIES ${SAGA_INCLUDE_DIR})
endif(SAGA_FOUND)

# muparser (optional for generic source spectra)
find_package(muParser)
if(MUPARSER_FOUND)
  list(APPEND CRPROPA_EXTRA_INCLUDES ${MUPARSER_INCLUDE_DIR})
  list(APPEND CRPROPA_EXTRA_LIBRARIES ${MUPARSER_LIBRARY})
  add_definitions (-DCRPROPA_HAVE_MUPARSER)
  list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_MUPARSER)
  list(APPEND SWIG_INCLUDE_DIRECTORIES ${MUPARSER_INCLUDE_DIR})
endif(MUPARSER_FOUND)

# zlib / zstream-cpp (optional for compressed output files)
find_package(ZLIB)
if(ZLIB_FOUND)
  list(APPEND CRPROPA_EXTRA_INCLUDES ${ZLIB_INCLUDE_DIRS})
  list(APPEND CRPROPA_EXTRA_INCLUDES "libs/zstream-cpp")
  list(APPEND CRPROPA_EXTRA_LIBRARIES ${ZLIB_LIBRARIES})
  add_definitions (-DCRPROPA_HAVE_ZLIB)
  list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_ZLIB)
  list(APPEND SWIG_INCLUDE_DIRECTORIES ${ZLIB_INCLUDE_DIRS})
endif(ZLIB_FOUND)

# HDF5 (optional for HDF5 output files)
option(ENABLE_HDF5 "HDF5 Support" ON)
if(ENABLE_HDF5)
  find_package( HDF5 COMPONENTS C )
  if(HDF5_FOUND)
    # parallel version requires MPI which is not supported in CRPropa
    if(NOT HDF5_IS_PARALLEL)
      list(APPEND CRPROPA_EXTRA_INCLUDES ${HDF5_INCLUDE_DIRS})
      list(APPEND CRPROPA_EXTRA_LIBRARIES ${HDF5_LIBRARIES})
      add_definitions (-DCRPROPA_HAVE_HDF5)
      list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_HDF5)
      list(APPEND SWIG_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
      #string(REPLACE " " " -I" HDF5_INCLUDE_DIRS_SWIG ${HDF5_INCLUDE_DIRS})
      #list(APPEND CRPROPA_SWIG_DEFINES -I${HDF5_INCLUDE_DIRS_SWIG})
    endif(NOT HDF5_IS_PARALLEL)
  endif(HDF5_FOUND)
endif(ENABLE_HDF5)


# ----------------------------------------------------------------------------
# Fix Apple RPATH
# ----------------------------------------------------------------------------
set(CMAKE_MACOSX_RPATH 1)
option(USE_ABSOLUTE_RPATH "Add absolute rpath to all libraries and executables" ON)
if(APPLE OR USE_ABSOLUTE_RPATH)
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(ABSOLUTE_RPATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
  if(NOT IS_ABSOLUTE ${ABSOLUTE_RPATH})
    set(ABSOLUTE_RPATH ${CMAKE_BINARY_DIR}/${ABSOLUTE_RPATH})
  endif()

  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${ABSOLUTE_RPATH}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    message(STATUS "Use absolute RPATH ${ABSOLUTE_RPATH}")
    set(CMAKE_INSTALL_RPATH "${ABSOLUTE_RPATH}")
  endif("${isSystemDir}" STREQUAL "-1")
endif()

# Fix for linker error on mac, Issue: #147
if(APPLE)
  add_definitions(-arch x86_64)
endif(APPLE)

# ----------------------------------------------------------------------------
# Download data files (interaction data, masses, decay data ...)
# ----------------------------------------------------------------------------
OPTION(DOWNLOAD_DATA "Download CRPropa data files" ON)
set(CRPROPA_DATAFILE_VER "2023-10-20")
if(DOWNLOAD_DATA)
  message("-- Downloading data files from sciebo ~ 73 MB")
  file(DOWNLOAD
  https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM
    ${CMAKE_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM 
    USERPWD "3juW9sntQX2IWBS")
  file(STRINGS ${CMAKE_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM DATA_CHECKSUM LIMIT_COUNT 1 LENGTH_MINIMUM 32 LENGTH_MAXIMUM 32)
  file(DOWNLOAD
  https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz
    ${CMAKE_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz
    USERPWD "3juW9sntQX2IWBS"
    EXPECTED_MD5 "${DATA_CHECKSUM}")
  message("-- Extracting data file")
else()
  message("-- Downloading of data file disabled")
endif(DOWNLOAD_DATA)
if(EXISTS ${CMAKE_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz)
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}/data ${CMAKE_BINARY_DIR}/data/ WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}/ WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
else()
  message(WARNING "CRPropa data file not found at ${CMAKE_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz
CRPropa should compile, but will likely not work properly! Please install data file manually, or use the automatic download which is enabled by default.")
endif()

# ----------------------------------------------------------------------------
# Library and Binary
# ----------------------------------------------------------------------------
file(GLOB_RECURSE CRPROPA_INCLUDES RELATIVE ${CMAKE_SOURCE_DIR} include/*.h)
include_directories(include ${CRPROPA_EXTRA_INCLUDES})

add_library(crpropa SHARED
  src/base64.cpp
  src/Candidate.cpp
  src/Clock.cpp
  src/Common.cpp
  src/Cosmology.cpp
  src/EmissionMap.cpp
  src/Geometry.cpp
  src/GridTools.cpp
  src/Module.cpp
  src/ModuleList.cpp
  src/ParticleID.cpp
  src/ParticleMass.cpp
  src/ParticleState.cpp
  src/PhotonBackground.cpp
  src/PhotonPropagation.cpp
  src/ProgressBar.cpp
  src/Random.cpp
  src/Source.cpp
  src/Variant.cpp
  src/module/AdiabaticCooling.cpp
  src/module/Acceleration.cpp
  src/module/Boundary.cpp
  src/module/BreakCondition.cpp
  src/module/CandidateSplitting.cpp
  src/module/DiffusionSDE.cpp
  src/module/EMCascade.cpp
  src/module/EMDoublePairProduction.cpp
  src/module/EMInverseComptonScattering.cpp
  src/module/EMPairProduction.cpp
  src/module/EMTripletPairProduction.cpp
  src/module/ElasticScattering.cpp
  src/module/ElectronPairProduction.cpp
  src/module/HDF5Output.cpp
  src/module/MomentumDiffusion.cpp
  src/module/NuclearDecay.cpp
  src/module/Observer.cpp
  src/module/Output.cpp
  src/module/OutputShell.cpp
  src/module/ParticleCollector.cpp
  src/module/PhotoDisintegration.cpp
  src/module/PhotoPionProduction.cpp
  src/module/PhotonEleCa.cpp
  src/module/PhotonOutput1D.cpp
  src/module/PropagationBP.cpp
  src/module/PropagationCK.cpp
  src/module/Redshift.cpp
  src/module/RestrictToRegion.cpp
  src/module/SimplePropagation.cpp
  src/module/SynchrotronRadiation.cpp
  src/module/TextOutput.cpp
  src/module/Tools.cpp
  src/magneticField/ArchimedeanSpiralField.cpp
  src/magneticField/JF12Field.cpp
  src/magneticField/JF12FieldSolenoidal.cpp
  src/magneticField/MagneticField.cpp
  src/magneticField/MagneticFieldGrid.cpp
  src/magneticField/PolarizedSingleModeMagneticField.cpp
  src/magneticField/PT11Field.cpp
  src/magneticField/turbulentField/GridTurbulence.cpp
  src/magneticField/turbulentField/HelicalGridTurbulence.cpp
  src/magneticField/turbulentField/PlaneWaveTurbulence.cpp
  src/magneticField/turbulentField/SimpleGridTurbulence.cpp
  src/magneticField/TF17Field.cpp
  src/magneticField/CMZField.cpp
  src/advectionField/AdvectionField.cpp
  src/massDistribution/ConstantDensity.cpp
  src/massDistribution/Cordes.cpp
  src/massDistribution/Ferriere.cpp
  src/massDistribution/Massdistribution.cpp
  src/massDistribution/Nakanishi.cpp

  ${CRPROPA_EXTRA_SOURCES}
)
target_link_libraries(crpropa ${CRPROPA_EXTRA_LIBRARIES})

#------------------------------------------------------------------
# Doxygen ; xml data is used for sphinx site and python docstrings
#------------------------------------------------------------------
option(BUILD_DOC "Enable building documentation from the source code" OFF)

if(BUILD_DOC)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    MESSAGE(STATUS "Found Doxygen to build documentation")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/DoxygenLayout.xml ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml COPYONLY)
    add_custom_target(doxy ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)

    set_source_files_properties(${CMAKE_BINARY_DIR}/docstrings_from_doxy.i PROPERTIES GENERATED true )
    add_custom_target(docstrings_from_doxy
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/libs/doxy2swig/doxy2swig.py ${CMAKE_CURRENT_BINARY_DIR}/xml/index.xml -o ${CMAKE_CURRENT_BINARY_DIR}/docstrings_from_doxy.i
      DEPENDS doxy
    )
  else()
    MESSAGE(STATUS "Doxygen not found: cannot build documentation")
  endif(DOXYGEN_FOUND)
else()
  MESSAGE(STATUS "Build of documentation disabeled. Enable with BUILD_DOC=On")
endif(BUILD_DOC)



# ----------------------------------------------------------------------------
# Python
# ----------------------------------------------------------------------------
option(ENABLE_PYTHON "Create python library via SWIG" ON)
find_package(PythonInterp)
find_package(PythonLibs)

if(ENABLE_PYTHON AND PYTHONLIBS_FOUND)

  find_package(SWIG 3.0 REQUIRED)
  include(python/Python.cmake)
  include_directories(${PYTHON_INCLUDE_PATH})

  if(SWIG_VERSION VERSION_GREATER 4.0)
    # Use swig 4 builtin doxygen instead of external program
    list(APPEND CRPROPA_SWIG_DEFINES -doxygen)
  else()
    if(BUILD_DOC AND DOXYGEN_FOUND)
      LIST(APPEND CRPROPA_SWIG_INPUTS ${CMAKE_BINARY_DIR}/docstrings_from_doxy.i)
      list(APPEND CRPROPA_SWIG_DEFINES -DWITH_DOXYGEN)
      list(APPEND SWIG_INCLUDE_DIRECTORIES ${CMAKE_BINARY_DIR})
    endif(BUILD_DOC AND DOXYGEN_FOUND)
  endif(SWIG_VERSION VERSION_GREATER 4.0)

  # use SWIG's builtin (maybe should be abandoned if everything will be tested and working with -builtin)
  option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" ON)

  if(ENABLE_SWIG_BUILTIN)
    set(BUILTIN "-builtin")
    set(PY3 "-py3")
    if(SWIG_VERSION VERSION_LESS 4.0.2)
        message(WARNING, "The SWIG builtin option should not be used with SWIG version below 4.0.2 due to https://github.com/swig/swig/issues/1595")
    endif()
  else(ENABLE_SWIG_BUILTIN)
    set(BUILTIN "")
    set(PY3 "")
  endif(ENABLE_SWIG_BUILTIN)

  if(PYTHON_VERSION_STRING VERSION_GREATER 3.0)
    list(APPEND CRPROPA_SWIG_DEFINES -DSWIG_PYTHON3)
  endif(PYTHON_VERSION_STRING VERSION_GREATER 3.0)

  # tries to import numpy
  execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/checkNumpy.py" OUTPUT_VARIABLE numpyIncludePath)
  if(numpyIncludePath)
    MESSAGE(STATUS "Found numpy headers in " ${numpyIncludePath})
    SET(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
    list(APPEND CRPROPA_SWIG_DEFINES -DWITHNUMPY)
    include_directories(${numpyIncludePath})
  else(numpyIncludePath)
    MESSAGE(STATUS "Numpy not found.")
  endif(numpyIncludePath)

  set(SWIG_INCLUDES)
  foreach(p in ${SWIG_INCLUDE_DIRECTORIES})
    list(APPEND SWIG_INCLUDES -I${p})
  endforeach()

  file(GLOB_RECURSE CRPROPA_SWIG_INPUTS python/*.i)
  set_source_files_properties( ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx PROPERTIES GENERATED true )
  add_custom_target(crpropa-swig-wrapper
    COMMAND swig ${BUILTIN} -c++ -python ${PY3} -I${CMAKE_SOURCE_DIR}/include -I${CMAKE_SOURCE_DIR}/libs/HepPID/include ${SWIG_INCLUDES} ${CRPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/python/crpropa${BUILTIN}.i
    DEPENDS ${CRPROPA_SWIG_INPUTS} ${CRPROPA_INCLUDES} )

  if(BUILD_DOC AND DOXYGEN_FOUND)
    add_dependencies(crpropa-swig-wrapper docstrings_from_doxy)
  endif(BUILD_DOC AND DOXYGEN_FOUND)

  add_library(crpropa-swig MODULE ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx)
  set_target_properties(crpropa-swig PROPERTIES PREFIX "")
  # disable warnings on automatically generated interface code
  set_target_properties(crpropa-swig PROPERTIES COMPILE_FLAGS "-w")
  set_target_properties(crpropa-swig PROPERTIES OUTPUT_NAME "_crpropa")
  target_link_libraries(crpropa-swig crpropa ${PYTHON_LIBRARIES})
  add_dependencies(crpropa-swig crpropa-swig-wrapper)

  install(DIRECTORY "${CMAKE_SOURCE_DIR}/python/crpropa" DESTINATION "${PYTHON_SITE_PACKAGES}")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/crpropa.py" DESTINATION "${PYTHON_SITE_PACKAGES}/crpropa")
  install(TARGETS crpropa-swig LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/crpropa")
  install(FILES ${CRPROPA_SWIG_INPUTS} DESTINATION share/crpropa/swig_interface)
endif(ENABLE_PYTHON AND PYTHONLIBS_FOUND)


# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
add_definitions(-DCRPROPA_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS crpropa DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/data/ DESTINATION share/crpropa/ PATTERN ".git" EXCLUDE)

install(DIRECTORY libs/kiss/include/ DESTINATION include)

# ------------------------------------------------------------------
# Documentation
# ------------------------------------------------------------------

find_program(SPHINX_EXECUTABLE
             NAMES sphinx-build
             DOC "Path to sphinx-build executable")

if(BUILD_DOC AND NOT SPHINX_EXECUTABLE)
  MESSAGE(STATUS "Sphinx not found: cannot build documentation")
endif(BUILD_DOC AND NOT SPHINX_EXECUTABLE)

if(BUILD_DOC AND DOXYGEN_FOUND AND SPHINX_EXECUTABLE)
  MESSAGE(STATUS "Found Sphinx to build documentation")
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Sphinx
                                    "Failed to find sphinx-build executable"
                                    SPHINX_EXECUTABLE)
  set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/doc)
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/doc)
  add_custom_target(doc
                    ${SPHINX_EXECUTABLE} -b html
                    -Dbreathe_projects.CRPropa3=${CMAKE_CURRENT_BINARY_DIR}/xml
                    ${SPHINX_SOURCE} ${SPHINX_BUILD}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating documentation with Sphinx")
  add_dependencies(doc doxy)
endif(BUILD_DOC AND DOXYGEN_FOUND AND SPHINX_EXECUTABLE)


# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------
if(ENABLE_TESTING)
  enable_testing()
  add_executable(testCore test/testCore.cpp)
  target_link_libraries(testCore  crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testCore testCore)

  add_executable(testVector3 test/testVector3.cpp)
  target_link_libraries(testVector3 crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testVector3 testVector3)

  add_executable(testModuleList test/testModuleList.cpp)
  target_link_libraries(testModuleList crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testModuleList testModuleList)

  add_executable(testMagneticField test/testMagneticField.cpp)
  target_link_libraries(testMagneticField crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testMagneticField testMagneticField)

  add_executable(testTurbulentField test/testTurbulentField.cpp)
  target_link_libraries(testTurbulentField crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testTurbulentField testTurbulentField)

  add_executable(testAdvectionField test/testAdvectionField.cpp)
  target_link_libraries(testAdvectionField crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testAdvectionField testAdvectionField)

  add_executable(testDensity test/testDensity.cpp)
  target_link_libraries(testDensity crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testDensity testDensity)

  add_executable(testDINT test/testDINT.cpp)
  target_link_libraries(testDINT crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testDINT testDINT)

  add_executable(testPropagation test/testPropagation.cpp)
  target_link_libraries(testPropagation crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testPropagation testPropagation)

  add_executable(testBreakCondition test/testBreakCondition.cpp)
  target_link_libraries(testBreakCondition crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testBreakCondition testBreakCondition)

  add_executable(testInteraction test/testInteraction.cpp)
  target_link_libraries(testInteraction crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testInteraction testInteraction)

  add_executable(testSource test/testSource.cpp)
  target_link_libraries(testSource crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testSource testSource)

  add_executable(testOutput test/testOutput.cpp)
  target_link_libraries(testOutput crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testOutput testOutput)

  add_executable(testFunctionalGroups test/testFunctionalGroups.cpp)
  target_link_libraries(testFunctionalGroups crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testFunctionalGroups testFunctionalGroups)

  add_executable(testAdiabaticCooling test/testAdiabaticCooling.cpp)
  target_link_libraries(testAdiabaticCooling crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testAdiabaticCooling testAdiabaticCooling)

  add_executable(testCandidateSplitting test/testCandidateSplitting.cpp)
  target_link_libraries(testCandidateSplitting crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
  add_test(testCandidateSplitting testCandidateSplitting)

  if(WITH_GALACTIC_LENSES)
    add_executable(testGalacticMagneticLens test/testMagneticLens.cpp)
    target_link_libraries(testGalacticMagneticLens crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
    add_test(testGalacticMagneticLens testGalacticMagneticLens)
  endif(WITH_GALACTIC_LENSES)

  # python tests
  if(ENABLE_PYTHON AND PYTHONLIBS_FOUND)
    CONFIGURE_FILE(test/testMagneticLensPythonInterface.py testMagneticLensPythonInterface.py COPYONLY)
    if(numpyIncludePath AND WITH_GALACTIC_LENSES)
      add_test(testMagneticLensPythonInterface ${PYTHON_EXECUTABLE} testMagneticLensPythonInterface.py)
    endif(numpyIncludePath AND WITH_GALACTIC_LENSES)

    CONFIGURE_FILE(test/testSimulationExecution.py testSimulationExecution.py COPYONLY)
    add_test(testSimulationExecution ${PYTHON_EXECUTABLE} testSimulationExecution.py)

    CONFIGURE_FILE(test/testDiffusionSDE.py testDiffusionSDE.py COPYONLY)
    add_test(testDiffusionSDE ${PYTHON_EXECUTABLE} testDiffusionSDE.py)

    CONFIGURE_FILE(test/testMomentumDiffusion.py testMomentumDiffusion.py COPYONLY)
    add_test(testMomentumDiffusion ${PYTHON_EXECUTABLE} testMomentumDiffusion.py)

    CONFIGURE_FILE(test/testPythonExtension.py testPythonExtension.py COPYONLY)
    add_test(testPythonExtension ${PYTHON_EXECUTABLE} testPythonExtension.py)
  endif(ENABLE_PYTHON AND PYTHONLIBS_FOUND)

endif(ENABLE_TESTING)
