cmake_minimum_required(VERSION 3.14)

project(CRPropa Fortran C CXX)
set(CRPROPA_RELEASE_VERSION 3.3.0+)  # Update for new release
set(CRPROPA_DATAFILE_VER "2026-07-08")  # Version of data files

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CRPROPA_EXTRA_SOURCES)
set(CRPROPA_EXTRA_INCLUDES)
set(CRPROPA_EXTRA_LIBRARIES)
set(CRPROPA_SWIG_DEFINES)
set(CRPROPA_SWIG_INPUTS)
set(CRPROPA_BUNDLED_LIBRARY_TARGETS)
set(CRPROPA_BUNDLED_INCLUDE_DIRS)

if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
  message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)

set(SIMD_EXTENSIONS "none" CACHE STRING "Choose which of the SIMD instruction set extensions your target CPU supports. Possible values are \"native\" (use everything that's supported by the CPU you're building on), \"none\", \"avx\", and \"avx+fma\".")

if(SIMD_EXTENSIONS STREQUAL "none")
  message("With SIMD_EXTENSIONS \"none\" tricubic interpolation of vector grids is not possible. You should set SIMD_EXTENSION to a compatible value (\"avx\", \"avx+fma\", or -- depending on the build CPU -- \"native\"). Trilinear interpolation is still possible.")
else()
  add_definitions(-DHAVE_SIMD)
endif()

if(SIMD_EXTENSIONS STREQUAL "avx")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -msse4.1 -msse4.2 -mavx" )
elseif(SIMD_EXTENSIONS STREQUAL "avx+fma")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -msse4.1 -msse4.2 -mavx -mfma" )
elseif(SIMD_EXTENSIONS STREQUAL "native")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -march=native" )
elseif(NOT SIMD_EXTENSIONS STREQUAL "none")
  message(SEND_ERROR "SIMD_EXTENSIONS must have one of these values: \"native\", \"none\", \"avx\", or \"avx+fma\".")
endif()

set(FAST_WAVES OFF CACHE BOOL "Enable SIMD optimizations for PlaneWaveTurbulence. Requires that SIMD_EXTENSIONS is set to a compatible value (\"avx\", \"avx+fma\", or -- depending on the build CPU -- \"native\"). The build will fail if this option is set but SIMD_EXTENSIONS is not compatible, so you can adjust your settings.")
if(FAST_WAVES)
  if(SIMD_EXTENSIONS STREQUAL "none")
    message(SEND_ERROR "Enabling FAST_WAVES while SIMD_EXTENSIONS is set to \"none\" won't work, since SIMD is required for FAST_WAVES. Please try disabling FAST_WAVES, or try using another option for SIMD_EXTENSIONS.")
  else()
    add_definitions(-DFAST_WAVES)
  endif()
endif(FAST_WAVES)

# Add build type for profiling
set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -ggdb -fno-omit-frame-pointer")

# Enable extra warnings on debug builds
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wextra")

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
option(ENABLE_GIT "Embedding information about CRPropa version from git" On)

if(ENABLE_GIT)
  include(GetGitRevisionDescription)
  get_git_head_revision(GIT_REFSPEC GIT_SHA1)
  git_describe(GIT_DESC --tags)
endif(ENABLE_GIT)
if(NOT ENABLE_GIT OR ("${GIT_SHA1}" STREQUAL "GITDIR-NOTFOUND"))
  set(GIT_REFSPEC "")
  set(GIT_SHA1 "")
  set(GIT_DESC "${CRPROPA_RELEASE_VERSION}-no-git")
endif()

message(STATUS "CRPropa version: ${GIT_DESC} ${GIT_SHA1} ${GIT_REFSPEC}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/Version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/src/Version.cpp" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/crpropa/Version.h.in" "${CMAKE_CURRENT_BINARY_DIR}/include/crpropa/Version.h" @ONLY)
list(APPEND CRPROPA_EXTRA_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/src/Version.cpp")
list(APPEND CRPROPA_EXTRA_INCLUDES "${CMAKE_CURRENT_BINARY_DIR}/include")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# googletest (provided, see https://code.google.com/p/googletest/wiki/FAQ
#     Why is it not recommended use a pre-compiled copy of Google Test?)
option(ENABLE_TESTING "Build tests and enable test target" On)
if(ENABLE_TESTING)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/gtest/include)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/libs/gtest)
  list(APPEND CRPROPA_BUNDLED_LIBRARY_TARGETS gtest gtest_main)
  list(APPEND CRPROPA_BUNDLED_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/libs/gtest/include)

  # temporary workaround for newer clang versions due to its handling of unwinding
  #   see: https://github.com/google/googletest/issues/3062
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
      add_definitions(-DCRPROPA_TESTS_SKIP_EXCEPTIONS)
    endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

  if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
  endif(APPLE)
endif(ENABLE_TESTING)

#
option(ENABLE_COVERAGE "Add code for coverage analysis." Off)
if(ENABLE_COVERAGE)
  find_program(LCOV_PATH lcov)
  find_program(GENHTML_PATH genhtml)
  if(LCOV_PATH AND GENHTML_PATH)
    message("Enabling coverage report via $make test && make coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    list(APPEND CRPROPA_EXTRA_LIBRARIES "-lgcov")
    list(APPEND CRPROPA_EXTRA_LIBRARIES "-fprofile-arcs")

    if(ENABLE_TESTING)
      set(COVERAGE_LIBS "-coverage -lgcov")
      add_custom_target(coverage_clean
        COMMAND ${LCOV_PATH} --directory ${CMAKE_CURRENT_BINARY_DIR} --zerocounters
      )
      add_custom_target(coverage
        # generate coverage data
        COMMAND ${LCOV_PATH} --directory ${CMAKE_CURRENT_BINARY_DIR} --capture --output-file ${CMAKE_CURRENT_BINARY_DIR}/coverage.info --ignore-errors negative,inconsistent,unused VERBATIM
        # clean external libs
        COMMAND ${LCOV_PATH} --remove ${CMAKE_CURRENT_BINARY_DIR}/coverage.info "/usr/include/*" "/usr/lib/*" "*/libs/gtest/*" "*/libs/zstream-cpp/*" "*/build/*" -o ${CMAKE_CURRENT_BINARY_DIR}/coverage.info.cleaned --ignore-errors negative,inconsistent,unused VERBATIM
        # Generate html output
        COMMAND ${GENHTML_PATH} -o ${CMAKE_CURRENT_BINARY_DIR}/coverageReport ${CMAKE_CURRENT_BINARY_DIR}/coverage.info.cleaned VERBATIM
        COMMAND echo "Generated coverage report in ${CMAKE_CURRENT_BINARY_DIR}/coverageReport/index.html"
      )
    endif(ENABLE_TESTING)

  else()
    if(NOT LCOV_PATH)
      message(WARNING "lcov not found, coverage report generation  not possible!")
    endif(NOT LCOV_PATH)
    if(NOT GENHTML_PATH)
      message(WARNING "genhtml not found, coverage report generation not possible!")
    endif(NOT GENHTML_PATH)
  endif(LCOV_PATH AND GENHTML_PATH)
endif(ENABLE_COVERAGE)

# kiss (provided)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/libs/kiss)
list(APPEND CRPROPA_EXTRA_LIBRARIES kiss)
list(APPEND CRPROPA_EXTRA_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/libs/kiss/include)
list(APPEND CRPROPA_BUNDLED_LIBRARY_TARGETS kiss)
list(APPEND CRPROPA_BUNDLED_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/libs/kiss/include)

# nanoflann (provided)
list(APPEND CRPROPA_EXTRA_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/libs/nanoflann/include)
list(APPEND CRPROPA_BUNDLED_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/libs/nanoflann/include)

# HepID (provided)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/libs/HepPID)
list(APPEND CRPROPA_EXTRA_LIBRARIES HepPID)
list(APPEND CRPROPA_EXTRA_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/libs/HepPID/include)
list(APPEND CRPROPA_BUNDLED_LIBRARY_TARGETS HepPID)
list(APPEND CRPROPA_BUNDLED_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/libs/HepPID/include)

# SOPHIA (provided)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/libs/sophia)
list(APPEND CRPROPA_EXTRA_LIBRARIES sophia gfortran)
list(APPEND CRPROPA_EXTRA_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/libs/sophia)

# OpenMP (optional for shared memory multiprocessing)
option(ENABLE_OPENMP "OpenMP for multithreading" On)
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
option(ENABLE_QUIMBY "Quimby Support" Off)
if(ENABLE_QUIMBY)
  find_package(Quimby)
  if(QUIMBY_FOUND)
    list(APPEND CRPROPA_EXTRA_INCLUDES ${QUIMBY_INCLUDE_DIR})
    list(APPEND CRPROPA_EXTRA_LIBRARIES ${QUIMBY_LIBRARY})
    add_definitions(-DCRPROPA_HAVE_QUIMBY)
    list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_QUIMBY)
    list(APPEND SWIG_INCLUDE_DIRECTORIES ${QUIMBY_INCLUDE_DIR}/../share/quimby)
    list(APPEND SWIG_INCLUDE_DIRECTORIES ${QUIMBY_INCLUDE_DIR})
  endif(QUIMBY_FOUND)
endif(ENABLE_QUIMBY)

# muparser (optional for generic source spectra)
find_package(muParser)
if(MUPARSER_FOUND)
  list(APPEND CRPROPA_EXTRA_INCLUDES ${MUPARSER_INCLUDE_DIR})
  list(APPEND CRPROPA_EXTRA_LIBRARIES ${MUPARSER_LIBRARY})
  add_definitions(-DCRPROPA_HAVE_MUPARSER)
  list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_MUPARSER)
  list(APPEND SWIG_INCLUDE_DIRECTORIES ${MUPARSER_INCLUDE_DIR})
endif(MUPARSER_FOUND)

# zlib / zstream-cpp (optional for compressed output files)
find_package(ZLIB)
if(ZLIB_FOUND)
	list(APPEND CRPROPA_BUNDLED_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/libs/zstream-cpp/include)
  list(APPEND CRPROPA_EXTRA_LIBRARIES ${ZLIB_LIBRARIES})
  list(APPEND CRPROPA_EXTRA_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/libs/zstream-cpp/include")
	add_definitions(-DCRPROPA_HAVE_ZLIB)
	list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_ZLIB)
	list(APPEND SWIG_INCLUDE_DIRECTORIES ${ZLIB_INCLUDE_DIRS})
endif(ZLIB_FOUND)

# HDF5 (optional for HDF5 output files)
option(ENABLE_HDF5 "HDF5 Support" On)
if(ENABLE_HDF5)
  find_package(HDF5 COMPONENTS C)
  if(HDF5_FOUND)
    # parallel version requires MPI which is not supported in CRPropa
    if(NOT HDF5_IS_PARALLEL)
      list(APPEND CRPROPA_EXTRA_INCLUDES ${HDF5_INCLUDE_DIRS})
      list(APPEND CRPROPA_EXTRA_LIBRARIES ${HDF5_LIBRARIES})
      add_definitions (-DCRPROPA_HAVE_HDF5)
      list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_HDF5)
      list(APPEND SWIG_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
    endif(NOT HDF5_IS_PARALLEL)
  endif(HDF5_FOUND)
endif(ENABLE_HDF5)


# ----------------------------------------------------------------------------
# Fix Apple RPATH
# ----------------------------------------------------------------------------

# for relative paths GNU-style
include(GNUInstallDirs)

file(RELATIVE_PATH _rpath_rel "${CMAKE_INSTALL_FULL_BINDIR}" "${CMAKE_INSTALL_FULL_LIBDIR}")
if(APPLE)
  set(_rpath_origin "@loader_path")
else()
  set(_rpath_origin "$ORIGIN")
endif()

# CMake emits the correct platform separator itself with lib beside the binary; second: bin -> lib.
set(CMAKE_INSTALL_RPATH "${_rpath_origin}" "${_rpath_origin}/${_rpath_rel}")

# Build tree keeps its own rpath so tests run against freshly built libs.
set(CMAKE_BUILD_WITH_INSTALL_RPATH False)




# ----------------------------------------------------------------------------
# Download data files (interaction data, masses, decay data ...)
# ----------------------------------------------------------------------------

option(DOWNLOAD_DATA "Download CRPropa data files" On)
if(DOWNLOAD_DATA)
  message("-- Downloading data files from sciebo ~ 260 MB")
  file(DOWNLOAD
  https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM
    ${CMAKE_CURRENT_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM
    USERPWD "Nag8KQLoCjpG5es")
  file(STRINGS ${CMAKE_CURRENT_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM DATA_CHECKSUM LIMIT_COUNT 1 LENGTH_MINIMUM 32 LENGTH_MAXIMUM 32)
  file(DOWNLOAD
    https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz
    ${CMAKE_CURRENT_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz
    USERPWD "Nag8KQLoCjpG5es"
    EXPECTED_MD5 "${DATA_CHECKSUM}")
else()
  message("-- Downloading of data file disabled")
endif(DOWNLOAD_DATA)
if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz)
  message("-- Extracting data file")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf data-${CRPROPA_DATAFILE_VER}.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory data-${CRPROPA_DATAFILE_VER}/data data/ WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory data-${CRPROPA_DATAFILE_VER}/ WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
else()
  message(WARNING "CRPropa data file not found at ${CMAKE_CURRENT_BINARY_DIR}/data-${CRPROPA_DATAFILE_VER}.tar.gz CRPropa should compile, but will likely not work properly! Please install data file manually, or use the automatic download which is enabled by default.")
endif()

# ----------------------------------------------------------------------------
# Library and Binary
# ----------------------------------------------------------------------------

file(GLOB_RECURSE CRPROPA_INCLUDES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} include/*.h)
# set CRPROPA_EXTRA_INCLUDES as cache variable to hand it over when used as subproject
set(CRPROPA_EXTRA_INCLUDES "${CRPROPA_EXTRA_INCLUDES}" CACHE STRING "The include paths of the extra libraries that are needed." FORCE)
mark_as_advanced(CRPROPA_EXTRA_INCLUDES)
include_directories(include ${CRPROPA_EXTRA_INCLUDES})

add_library(crpropa SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/src/base64.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Candidate.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Common.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Cosmology.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/EmissionMap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/GridTools.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/InteractionRates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Module.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleID.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleMass.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleState.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PhotonBackground.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ProgressBar.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Source.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Variant.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/AdiabaticCooling.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/Acceleration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/Boundary.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/BreakCondition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/CandidateSplitting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/DiffusionSDE.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/EMDoublePairProduction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/EMInverseComptonScattering.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/EMPairProduction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/EMTripletPairProduction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/ElasticScattering.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/ElectronPairProduction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/HDF5Output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/MomentumDiffusion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/NuclearDecay.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/Observer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/Output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/OutputShell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/ParticleCollector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/PhotoDisintegration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/PhotoPionProduction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/PhotonOutput1D.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/PropagationBP.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/PropagationCK.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/Redshift.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/RestrictToRegion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/SimplePropagation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/SynchrotronRadiation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/TextOutput.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/module/Tools.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/ArchimedeanSpiralField.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/JF12Field.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/JF12FieldSolenoidal.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/MagneticField.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/MagneticFieldGrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/PolarizedSingleModeMagneticField.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/PT11Field.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/turbulentField/GridTurbulence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/turbulentField/HelicalGridTurbulence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/turbulentField/PlaneWaveTurbulence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/turbulentField/SimpleGridTurbulence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/TF17Field.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/UF23Field.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/KST24Field.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/magneticField/CMZField.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/advectionField/AdvectionField.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/advectionField/TimeDependentAdvectionField.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/massDistribution/ConstantDensity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/massDistribution/Cordes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/massDistribution/Ferriere.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/massDistribution/Massdistribution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/massDistribution/Nakanishi.cpp

  ${CRPROPA_EXTRA_SOURCES}
)
target_link_libraries(crpropa ${CRPROPA_EXTRA_LIBRARIES})

#------------------------------------------------------------------
# Doxygen ; xml data is used for sphinx site and python docstrings
#------------------------------------------------------------------
option(BUILD_DOC "Enable building documentation from the source code" Off)

if(BUILD_DOC)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    message(STATUS "Found Doxygen to build documentation")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/DoxygenLayout.xml ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml COPYONLY)
    add_custom_target(doxy ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)

    set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/docstrings_from_doxy.i PROPERTIES GENERATED true )
    add_custom_target(docstrings_from_doxy
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/libs/doxy2swig/doxy2swig.py ${CMAKE_CURRENT_BINARY_DIR}/xml/index.xml -o ${CMAKE_CURRENT_BINARY_DIR}/docstrings_from_doxy.i
      DEPENDS doxy
    )
  else()
    message(STATUS "Doxygen not found: cannot build documentation")
  endif(DOXYGEN_FOUND)
else()
  message(STATUS "Build of documentation disabeled. Enable with BUILD_DOC=On")
endif(BUILD_DOC)

# Only add stdc++fs for older GCC (needed for Ubuntu 20)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    target_link_libraries(crpropa stdc++fs)
endif()

# ----------------------------------------------------------------------------
# Python
# ----------------------------------------------------------------------------
option(ENABLE_PYTHON "Create python library via SWIG" On)

set(Python_FIND_VIRTUALENV "FIRST")
set(Python_FIND_STRATEGY "LOCATION")
find_package(Python 3.0 COMPONENTS Interpreter Development NumPy)

if(ENABLE_PYTHON AND Python_FOUND)
  find_package(SWIG 3.0 REQUIRED)
  include(UseSWIG)
  include_directories(${Python_INCLUDE_DIRS})
  add_definitions(-DCRPROPA_HAVE_PYTHON)
  list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_PYTHON)

  # print Python info in detail
  message(STATUS "Python: Found!")
  message(STATUS "  version ${Python_VERSION}")
  message(STATUS "  executable: ${Python_EXECUTABLE}")
  message(STATUS "  libraries: ${Python_LIBRARIES}")
  message(STATUS "  headers: ${Python_INCLUDE_DIRS}")
  message(STATUS "  site packages: ${Python_SITELIB}")
  if(Python_Development_FOUND)
    message(STATUS "  development libraries: Found!")
  elseif(Python_Development_FOUND)
    message(STATUS "  development libraries: NOT found!")
  endif(Python_Development_FOUND)

  
  #  use Python_INSTALL_PACKAGE_DIR if provided; otherwise, install in Python_SITELIB
  set(Python_INSTALL_PACKAGE_DIR "${Python_SITELIB}" CACHE PATH "folder in which the python package is installed")
  message(STATUS "  package install directory: ${Python_INSTALL_PACKAGE_DIR}")
  
  # look for NumPy
  if(Python_NumPy_FOUND)
    set(CMAKE_SWIG_FLAGS -DWITHNUMPY)
    list(APPEND CRPROPA_SWIG_DEFINES -DWITHNUMPY)
    include_directories(${Python_NumPy_INCLUDE_DIRS})
    message(STATUS "NumPy: Found!")
    message(STATUS "  headers: ${Python_NumPy_INCLUDE_DIRS} (version ${Python_NumPy_VERSION})")
  elseif(Python_NumPy_FOUND)
    message(STATUS "NumPy: NOT found!")
  endif(Python_NumPy_FOUND)
  
  if(SWIG_VERSION VERSION_GREATER 4.0)
    # Use swig 4 builtin doxygen instead of external program
    list(APPEND CRPROPA_SWIG_DEFINES -doxygen)
  else()
    if(BUILD_DOC AND DOXYGEN_FOUND)
      list(APPEND CRPROPA_SWIG_INPUTS ${CMAKE_CURRENT_BINARY_DIR}/docstrings_from_doxy.i)
      list(APPEND CRPROPA_SWIG_DEFINES -DWITH_DOXYGEN)
      list(APPEND SWIG_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR})
    endif(BUILD_DOC AND DOXYGEN_FOUND)
  endif(SWIG_VERSION VERSION_GREATER 4.0)

  # use SWIG's builtin (maybe should be abandoned if everything will be tested and working with -builtin)
  option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" On)

  if(ENABLE_SWIG_BUILTIN)
    set(BUILTIN "-builtin")
    list(APPEND CRPROPA_SWIG_DEFINES "-builtin")
    if(SWIG_VERSION VERSION_LESS 4.0.2)
        message(WARNING, "The SWIG builtin option should not be used with SWIG version below 4.0.2 due to https://github.com/swig/swig/issues/1595")
    endif()
  else()
    set(BUILTIN "")
  endif(ENABLE_SWIG_BUILTIN)

  set(SWIG_INCLUDES)
  foreach(p in ${SWIG_INCLUDE_DIRECTORIES})
    list(APPEND SWIG_INCLUDES ${p})
  endforeach()

  file(GLOB_RECURSE CRPROPA_SWIG_INPUTS python/*.i)

  set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/python/crpropa${BUILTIN}.i
    PROPERTIES
    CPLUSPLUS ON
    SWIG_MODULE_NAME "crpropa"
    SWIG_FLAGS "${CRPROPA_SWIG_DEFINES};-dirprot;${CMAKE_SWIG_FLAGS}"
    DEPENDS "${CRPROPA_SWIG_INPUTS};${CRPROPA_INCLUDES}"
    OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
  
  swig_add_library(crpropa-swig TYPE SHARED LANGUAGE python
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/python/crpropa${BUILTIN}.i
  )
  
  target_include_directories(crpropa-swig PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/HepPID/include
    ${SWIG_INCLUDES}
  )
  
  target_link_libraries(crpropa-swig crpropa ${Python_LIBRARIES} ${Python_LIBRARY})
  set_target_properties(crpropa-swig PROPERTIES
    PREFIX ""
    OUTPUT_NAME "_crpropa"
    SUFFIX ".so"
  )

  if(APPLE)
    set_target_properties(crpropa-swig PROPERTIES
      BUILD_RPATH "$<TARGET_FILE_DIR:crpropa>"
      INSTALL_RPATH "@loader_path;${CMAKE_INSTALL_FULL_LIBDIR}"
    )
    set_property(TARGET crpropa-swig APPEND PROPERTY LINK_FLAGS "-flat_namespace -undefined suppress")
  endif(APPLE)
  

  if(BUILD_DOC AND DOXYGEN_FOUND)
    add_dependencies(crpropa-swig docstrings_from_doxy)
  endif(BUILD_DOC AND DOXYGEN_FOUND)

  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/crpropa" DESTINATION "${Python_INSTALL_PACKAGE_DIR}")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/crpropa.py" DESTINATION "${Python_INSTALL_PACKAGE_DIR}/crpropa")
  install(TARGETS crpropa-swig LIBRARY DESTINATION "${Python_INSTALL_PACKAGE_DIR}/crpropa")
  install(FILES ${CRPROPA_SWIG_INPUTS} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/crpropa/swig_interface)

endif(ENABLE_PYTHON AND Python_FOUND)

# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------

# makes the binary directory visible to common.h, needed to have a last resort when searching for data
add_definitions(-DCRPROPA_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
add_definitions(-DCRPROPA_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS crpropa ${CRPROPA_BUNDLED_LIBRARY_TARGETS}
	RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
	ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h")
if(DOWNLOAD_DATA)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data/ DESTINATION ${CMAKE_INSTALL_PREFIX}/share/crpropa/ PATTERN ".git" EXCLUDE)
endif()
foreach(CRPROPA_BUNDLED_INCLUDE_DIR IN LISTS CRPROPA_BUNDLED_INCLUDE_DIRS)
	install(DIRECTORY ${CRPROPA_BUNDLED_INCLUDE_DIR}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h" PATTERN "*.hh" PATTERN "*.hpp" PATTERN "*.hxx")
endforeach()

option(ENABLE_PYBIND11_STUBGEN "Generate Pybind11 stubs" On)
if(ENABLE_PYBIND11_STUBGEN)
  if(ENABLE_PYTHON AND Python_FOUND)
    find_program(PYBIND11_STUBGEN_EXECUTABLE NAMES "pybind11-stubgen")
    if(PYBIND11_STUBGEN_EXECUTABLE)
      message(STATUS "Found pybind11-stubgen at ${PYBIND11_STUBGEN_EXECUTABLE}")
      install(CODE "execute_process(COMMAND ${PYBIND11_STUBGEN_EXECUTABLE} --ignore-all-errors -o ${Python_SITELIB} crpropa)")
      install(CODE "execute_process(COMMAND echo -- Generated python stubs over pybind11-stubgen for crpropa)")
    endif(PYBIND11_STUBGEN_EXECUTABLE)
  endif(ENABLE_PYTHON AND Python_FOUND)
endif(ENABLE_PYBIND11_STUBGEN)



# ------------------------------------------------------------------
# Documentation
# ------------------------------------------------------------------

find_program(SPHINX_EXECUTABLE NAMES sphinx-build DOC "Path to sphinx-build executable")

if(BUILD_DOC AND NOT SPHINX_EXECUTABLE)
  message(STATUS "Sphinx not found: cannot build documentation")
endif(BUILD_DOC AND NOT SPHINX_EXECUTABLE)

if(BUILD_DOC AND DOXYGEN_FOUND AND SPHINX_EXECUTABLE)
  message(STATUS "Found Sphinx to build documentation")
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Sphinx "Failed to find sphinx-build executable" SPHINX_EXECUTABLE)
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/doc "${CMAKE_CURRENT_BINARY_DIR}/doc-source")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/xml "${CMAKE_CURRENT_BINARY_DIR}/doc-source/xml")
  set(SPHINX_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/doc-source")
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/doc)
  add_custom_target(doc
    ${SPHINX_EXECUTABLE} -b html
    -Dbreathe_projects.CRPropa3=xml
    ${SPHINX_SOURCE} ${SPHINX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Sphinx"
  )
  add_dependencies(doc doxy)
endif(BUILD_DOC AND DOXYGEN_FOUND AND SPHINX_EXECUTABLE)


# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------

if(ENABLE_TESTING)

  enable_testing()
  set(CRPROPA_TESTS_PATH ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "")  
  mark_as_advanced(CRPROPA_TESTS_PATH)

  function(addTest testName)
    add_executable(${testName} ${CMAKE_CURRENT_SOURCE_DIR}/test/${testName}.cpp)
    target_link_libraries(${testName} PRIVATE crpropa gtest gtest_main pthread ${COVERAGE_LIBS})
    add_test(NAME ${testName} COMMAND ${CRPROPA_TESTS_PATH}/${testName})
  endfunction()

  addTest(testCore)
  addTest(testVector3)
  addTest(testModuleList)
  addTest(testMagneticField)
  addTest(testTurbulentField)
  addTest(testAdvectionField)
  addTest(testDensity)
  addTest(testPropagation)
  addTest(testBreakCondition)
  addTest(testInteraction)
  addTest(testSource)
  addTest(testOutput)
  addTest(testFunctionalGroups)
  addTest(testAdiabaticCooling)
  addTest(testCandidateSplitting)

  # python tests
  if(ENABLE_PYTHON AND Python_FOUND)

    function(addTestPy testName)
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/${testName}.py ${CMAKE_CURRENT_BINARY_DIR}/${testName}.py COPYONLY)
      add_test(NAME ${testName} COMMAND ${Python_EXECUTABLE} ${CRPROPA_TESTS_PATH}/${testName}.py)
    endfunction()

    addTestPy(testSimulationExecution)
    addTestPy(testDiffusionSDE)
    addTestPy(testMomentumDiffusion)
    addTestPy(testPythonExtension)

  endif(ENABLE_PYTHON AND Python_FOUND)

endif(ENABLE_TESTING)
