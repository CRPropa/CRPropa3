if(APPLE)
	# rpath specific patches
	cmake_minimum_required(VERSION 2.8.12)
else(APPLE)
	# require > 2.8.8 for FILE DOWNLOAD fixes
	# allow < 2.8.12 for debian backports
	cmake_minimum_required(VERSION 2.8.11)
endif(APPLE)
project(CRPropa Fortran C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CRPROPA_EXTRA_SOURCES)
set(CRPROPA_EXTRA_INCLUDES)
set(CRPROPA_EXTRA_LIBRARIES)
set(CRPROPA_SWIG_DEFINES)

if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
	message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)

# Add build type for profiling
SET(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -ggdb -fno-omit-frame-pointer")

# Set default build-type to release to enable performance improvements
if (NOT CMAKE_BUILD_TYPE)
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
else(ENABLE_GIT)
	set(GIT_REFSPEC "")
	set(GIT_SHA1 "")
	set(GIT_DESC "3.0-no-git")
endif(ENABLE_GIT)

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

# kiss (provided)
add_subdirectory(libs/kiss)
list(APPEND CRPROPA_EXTRA_LIBRARIES kiss)
list(APPEND CRPROPA_EXTRA_INCLUDES libs/kiss/include)

# pugixml (provided)
add_subdirectory(libs/pugixml)
list(APPEND CRPROPA_EXTRA_LIBRARIES pugixml)
list(APPEND CRPROPA_EXTRA_INCLUDES libs/pugixml)

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
if (ENABLE_GALACTICMAGETICLENS)
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
endif (ENABLE_GALACTICMAGETICLENS)

# ROOT (optional for ROOT output)
option(ENABLE_ROOT "ROOT Output" ON)
if(ENABLE_ROOT)
	find_package(ROOT)
	if(ROOT_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CFLAGS}")
		set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ROOT_LIBS}")
		set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${ROOT_LIBS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ROOT_LIBS}")
		add_definitions(-DCRPROPA_HAVE_ROOT)
		list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_ROOT)
		list(APPEND CRPROPA_EXTRA_LIBRARIES Core RIO Tree)
	endif(ROOT_FOUND)
endif(ENABLE_ROOT)

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
if (ENABLE_QUIMBY)
	find_package(Quimby)
	if(QUIMBY_FOUND)
		list(APPEND CRPROPA_EXTRA_INCLUDES ${QUIMBY_INCLUDE_DIR})
		list(APPEND CRPROPA_EXTRA_LIBRARIES ${QUIMBY_LIBRARY})
		add_definitions (-DCRPROPA_HAVE_QUIMBY)
		list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_QUIMBY)
		list(APPEND CRPROPA_SWIG_DEFINES -I${QUIMBY_INCLUDE_DIR}/../share/quimby)
		list(APPEND CRPROPA_SWIG_DEFINES -I${QUIMBY_INCLUDE_DIR})
	endif(QUIMBY_FOUND)
endif(ENABLE_QUIMBY)

# SAGA (optinal for AMR grids)
find_package(SAGA)
if(SAGA_FOUND)
	list(APPEND CRPROPA_EXTRA_INCLUDES ${SAGA_INCLUDE_DIR})
	list(APPEND CRPROPA_EXTRA_LIBRARIES ${SAGA_LIBRARY})
	add_definitions (-DCRPROPA_HAVE_SAGA)
	list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_SAGA)
	list(APPEND CRPROPA_SWIG_DEFINES -I${SAGA_INCLUDE_DIR}/../python)
	list(APPEND CRPROPA_SWIG_DEFINES -I${SAGA_INCLUDE_DIR})
endif(SAGA_FOUND)

# muparser (optional for generic source spectra)
find_package(muParser)
if(MUPARSER_FOUND)
	list(APPEND CRPROPA_EXTRA_INCLUDES ${MUPARSER_INCLUDE_DIR})
	list(APPEND CRPROPA_EXTRA_LIBRARIES ${MUPARSER_LIBRARY})
	add_definitions (-DCRPROPA_HAVE_MUPARSER)
	list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_MUPARSER)
	list(APPEND CRPROPA_SWIG_DEFINES -I${MUPARSER_INCLUDE_DIR})
endif(MUPARSER_FOUND)

# zlib / zstream-cpp (optional for compressed output files)
find_package(ZLIB)
if(ZLIB_FOUND)
	list(APPEND CRPROPA_EXTRA_INCLUDES ${ZLIB_INCLUDE_DIRS})
	list(APPEND CRPROPA_EXTRA_INCLUDES "libs/zstream-cpp")
	list(APPEND CRPROPA_EXTRA_LIBRARIES ${ZLIB_LIBRARIES})
	add_definitions (-DCRPROPA_HAVE_ZLIB)
	list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_ZLIB)
	list(APPEND CRPROPA_SWIG_DEFINES -I${ZLIB_INCLUDE_DIRS})
endif(ZLIB_FOUND)

# HDF5 (optional for HDF5 output files)
option(ENABLE_HDF5 "HDF5 Support" ON)
if(ENABLE_HDF5)
	find_package( HDF5 COMPONENTS C )
	if(HDF5_FOUND)
		list(APPEND CRPROPA_EXTRA_INCLUDES ${HDF5_INCLUDE_DIRS})
		list(APPEND CRPROPA_EXTRA_LIBRARIES ${HDF5_LIBRARIES})
		add_definitions (-DCRPROPA_HAVE_HDF5)
		list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_HDF5)
		list(APPEND CRPROPA_SWIG_DEFINES -I${HDF5_INCLUDE_DIRS})
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


# ----------------------------------------------------------------------------
# Download data files (interaction data, masses, decay data ...)
# ----------------------------------------------------------------------------
message("-- Downloading data file from crpropa.desy.de ~ 35 MB")
file(DOWNLOAD
	https://crpropa.desy.de/images/8/83/Data_2016_05_23.tar.gz
	${CMAKE_BINARY_DIR}/data.tar.gz
	EXPECTED_MD5 6dd16796ba7517af562ab65a69fecac3)
message("-- Extracting data file")
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_BINARY_DIR}/data.tar.gz WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/data)


# ----------------------------------------------------------------------------
# Library and Binary
# ----------------------------------------------------------------------------
include_directories(include ${CRPROPA_EXTRA_INCLUDES})

add_library(crpropa SHARED
	src/Random.cpp
	src/Clock.cpp
	src/ModuleList.cpp
	src/Module.cpp
	src/Candidate.cpp
	src/EmissionMap.cpp
	src/ParticleState.cpp
	src/ParticleID.cpp
	src/ParticleMass.cpp
	src/ProgressBar.cpp
	src/Cosmology.cpp
	src/Source.cpp
	src/Common.cpp
	src/PhotonBackground.cpp
	src/PhotonPropagation.cpp
	src/GridTools.cpp
	src/XmlExecute.cpp
	src/module/BreakCondition.cpp
	src/module/Boundary.cpp
	src/module/Observer.cpp
	src/module/SimplePropagation.cpp
	src/module/PropagationCK.cpp
	src/module/ElectronPairProduction.cpp
	src/module/NuclearDecay.cpp
	src/module/PhotoPionProduction.cpp
	src/module/PhotoDisintegration.cpp
	src/module/Redshift.cpp
	src/module/Output.cpp
	src/module/TextOutput.cpp
	src/module/HDF5Output.cpp
	src/module/OutputShell.cpp
	src/module/OutputROOT.cpp
	src/module/OutputCRPropa2.cpp
	src/module/PhotonOutput1D.cpp
	src/module/PhotonEleCa.cpp
	src/module/EMCascade.cpp
	src/module/EMPairProduction.cpp
	src/module/EMDoublePairProduction.cpp
	src/module/EMTripletPairProduction.cpp
	src/module/EMInverseComptonScattering.cpp
	src/module/SynchrotronRadiation.cpp
	src/module/Tools.cpp
	src/magneticField/MagneticField.cpp
	src/magneticField/MagneticFieldGrid.cpp
	src/magneticField/JF12Field.cpp
	src/magneticField/PshirkovField.cpp
	${CRPROPA_EXTRA_SOURCES}
)
target_link_libraries(crpropa ${CRPROPA_EXTRA_LIBRARIES})

add_executable(crpropa-xmlrun src/main.cpp)
target_link_libraries(crpropa-xmlrun crpropa)


# ----------------------------------------------------------------------------
# Python
# ----------------------------------------------------------------------------
option(ENABLE_PYTHON "Create python library via SWIG" ON)
if(ENABLE_PYTHON)
	FIND_PACKAGE(SWIG REQUIRED)
	include(python/Python.cmake)
	include_directories(${PYTHON_INCLUDE_PATH})

	# use SWIG's builtin (maybe should be abandoned if everything will be tested and working with -builtin)
	option(ENABLE_SWIG_BUILTIN "Use SWIG's -builtin option" OFF)
	if(ENABLE_SWIG_BUILTIN)
		set(BUILTIN "-builtin")
	else(ENABLE_SWIG_BUILTIN)
		set(BUILTIN "")
	endif(ENABLE_SWIG_BUILTIN)

	# tries to import numpy
	execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/checkNumpy.py" OUTPUT_VARIABLE numpyIncludePath)
	if(numpyIncludePath)
	  MESSAGE(STATUS "Found numpy headers in " ${numpyIncludePath})
		SET(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
		list(APPEND CRPROPA_SWIG_DEFINES -DWITHNUMPY)
	  include_directories(${numpyIncludePath})
	else()
	  MESSAGE(STATUS "Numpy not found.")
	endif()

	file(GLOB_RECURSE CRPROPA_INCLUDES include/*.h)
	file(GLOB_RECURSE CRPROPA_SWIG_INPUTS python/*.i)
	set_source_files_properties( ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx PROPERTIES GENERATED true )
	add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx
		COMMAND swig ${BUILTIN} -c++ -python -I${CMAKE_SOURCE_DIR}/include ${CRPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/python/crpropa${BUILTIN}.i
				DEPENDS ${CRPROPA_SWIG_INPUTS} ${CRPROPA_INCLUDES} )

	add_library(crpropa-swig MODULE ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx)
	set_target_properties(crpropa-swig PROPERTIES PREFIX "")
	set_target_properties(crpropa-swig PROPERTIES OUTPUT_NAME "_crpropa")
	target_link_libraries(crpropa-swig crpropa ${PYTHON_LIBRARIES})

	install(DIRECTORY "${CMAKE_SOURCE_DIR}/python/crpropa" DESTINATION "${PYTHON_SITE_PACKAGES}")
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/crpropa.py" DESTINATION "${PYTHON_SITE_PACKAGES}/crpropa")
	install(TARGETS crpropa-swig LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/crpropa")
endif(ENABLE_PYTHON)


# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
add_definitions(-DCRPROPA_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS crpropa-xmlrun RUNTIME DESTINATION bin)
install(TARGETS crpropa DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/data/ DESTINATION share/crpropa/ PATTERN ".git" EXCLUDE)


# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------
if(ENABLE_TESTING)
	enable_testing()
	add_executable(testCore test/testCore.cpp)
	target_link_libraries(testCore crpropa gtest gtest_main pthread)
	add_test(testCore testCore)

	add_executable(testVector3 test/testVector3.cpp)
	target_link_libraries(testVector3 crpropa gtest gtest_main pthread)
	add_test(testVector3 testVector3)

	add_executable(testModuleList test/testModuleList.cpp)
	target_link_libraries(testModuleList crpropa gtest gtest_main pthread)
	add_test(testModuleList testModuleList)

	add_executable(testMagneticField test/testMagneticField.cpp)
	target_link_libraries(testMagneticField crpropa gtest gtest_main pthread)
	add_test(testMagneticField testMagneticField)

	add_executable(testDINT test/testDINT.cpp)
	target_link_libraries(testDINT crpropa gtest gtest_main pthread)
	add_test(testDINT testDINT)

	add_executable(testPropagation test/testPropagation.cpp)
	target_link_libraries(testPropagation crpropa gtest gtest_main pthread)
	add_test(testPropagation testPropagation)

	add_executable(testBreakCondition test/testBreakCondition.cpp)
	target_link_libraries(testBreakCondition crpropa gtest gtest_main pthread)
	add_test(testBreakCondition testBreakCondition)

	add_executable(testInteraction test/testInteraction.cpp)
	target_link_libraries(testInteraction crpropa gtest gtest_main pthread)
	add_test(testInteraction testInteraction)

	add_executable(testSource test/testSource.cpp)
	target_link_libraries(testSource crpropa gtest gtest_main pthread)
	add_test(testSource testSource)

	if(WITH_GALACTIC_LENSES)
		add_executable(testGalacticMagneticLens test/testMagneticLens.cpp)
		target_link_libraries(testGalacticMagneticLens crpropa gtest gtest_main pthread)
		add_test(testGalacticMagneticLens testGalacticMagneticLens)
	endif(WITH_GALACTIC_LENSES)

	if(ENABLE_PYTHON)
		CONFIGURE_FILE(test/testMagneticLensPythonInterface.py.in testMagneticLensPythonInterface.py)
		add_test(testMagneticLensPythonInterface ${PYTHON_EXECUTABLE} testMagneticLensPythonInterface.py)

		CONFIGURE_FILE(test/testSimulationExecution.py.in testSimulationExecution.py)
		add_test(testSimulationExecution ${PYTHON_EXECUTABLE} testSimulationExecution.py)

		# Pure python tests
		add_custom_command(TARGET crpropa POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/test/python/*.py ${CMAKE_CURRENT_BINARY_DIR}/)

		add_test(testPythonBoundaries ${PYTHON_EXECUTABLE} testBoundaries.py)
		add_test(testPythonCosmology ${PYTHON_EXECUTABLE} testCosmology.py)
		add_test(testPythonElectronPairProduction ${PYTHON_EXECUTABLE} testElectronPairProduction.py)
		add_test(testPythonEnergyLossLength ${PYTHON_EXECUTABLE} testEnergyLossLength.py)
		if(WITH_GALACTIC_LENSES)
			add_test(testPythonJF12Field ${PYTHON_EXECUTABLE} testJF12Field.py)
		endif(WITH_GALACTIC_LENSES)
		add_test(testPythonNuclearDecay ${PYTHON_EXECUTABLE} testNuclearDecay.py)
		add_test(testPythonNuclearMass ${PYTHON_EXECUTABLE} testNuclearMass.py)
		if(ENABLE_HDF5)
			add_test(testPythonOutputHDF5 ${PYTHON_EXECUTABLE} testOutputHDF5.py)
		endif(ENABLE_HDF5)
		add_test(testPythonPhotoPionProduction_pnRatio ${PYTHON_EXECUTABLE} testPhotoPionProduction_pnRatio.py)
		add_test(testPythonPropagationCK ${PYTHON_EXECUTABLE} testPropagationCK.py)
		add_test(testPythonPshirkovField ${PYTHON_EXECUTABLE} testPshirkovField.py)
		add_test(testPythonSourceComposition ${PYTHON_EXECUTABLE} testSourceComposition.py)
		add_test(testPythonTools ${PYTHON_EXECUTABLE} testTools.py)
		#add_test(testPythonTrackingObserver ${PYTHON_EXECUTABLE} testTrackingObserver.py) # lasts for too long? 10+min
		if(FFTW3F_FOUND)
			add_test(testPythonTurbulentDeflection ${PYTHON_EXECUTABLE} testTurbulentDeflection.py)
			add_test(testPythonTurbulentGrid ${PYTHON_EXECUTABLE} testTurbulentGrid.py)
		endif(FFTW3F_FOUND)
	endif(ENABLE_PYTHON)
endif(ENABLE_TESTING)
