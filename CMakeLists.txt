cmake_minimum_required(VERSION 2.6)
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

# Set default build-type to release to enable performance improvements 
if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
MESSAGE(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

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
SET(WITH_GALACTIC_LENSES FALSE)
if (ENABLE_GALACTICMAGETICLENS)
	FIND_PACKAGE(Boost)
	IF(Boost_FOUND)
		SET(WITH_GALACTIC_LENSES TRUE)

		list(APPEND CRPROPA_SWIG_DEFINES -DWITH_GALACTIC_LENSES)

		list(APPEND CRPROPA_EXTRA_SOURCES src/magneticLens/MagneticLens.cpp)
		list(APPEND CRPROPA_EXTRA_SOURCES src/magneticLens/ModelMatrix.cpp)
		list(APPEND CRPROPA_EXTRA_SOURCES	src/magneticLens/Pixelization.cpp)
		list(APPEND CRPROPA_EXTRA_SOURCES	src/magneticLens/ParticleMapsContainer.cpp)

		list(APPEND CRPROPA_EXTRA_INCLUDES ${Boost_INCLUDE_DIRS})
		
	ELSE(Boost_FOUND)
		MESSAGE (STATUS "Cannot find BOOST which is required for GalacticMagneticLenses.")
	ENDIF(Boost_FOUND)
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
find_package(Quimby)
if(QUIMBY_FOUND)
	list(APPEND CRPROPA_EXTRA_INCLUDES ${QUIMBY_INCLUDE_DIR})
	list(APPEND CRPROPA_EXTRA_LIBRARIES ${QUIMBY_LIBRARY})
	add_definitions (-DCRPROPA_HAVE_QUIMBY)
	list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_QUIMBY)
	list(APPEND CRPROPA_SWIG_DEFINES -I${QUIMBY_INCLUDE_DIR}/../share/quimby)
	list(APPEND CRPROPA_SWIG_DEFINES -I${QUIMBY_INCLUDE_DIR})
endif(QUIMBY_FOUND)

# muparser (optional for genercic source spectra)
find_package(muParser)
if(MUPARSER_FOUND)
	list(APPEND CRPROPA_EXTRA_INCLUDES ${MUPARSER_INCLUDE_DIR})
	list(APPEND CRPROPA_EXTRA_LIBRARIES ${MUPARSER_LIBRARY})
	add_definitions (-DCRPROPA_HAVE_MUPARSER)
	list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_MUPARSER)
	list(APPEND CRPROPA_SWIG_DEFINES -I${MUPARSER_INCLUDE_DIR})
endif(MUPARSER_FOUND)

# ----------------------------------------------------------------------------
# Downloading data files (interaction data, masses, decay data)
# ----------------------------------------------------------------------------
set(DATA_URL "http://crpropa.desy.de/images/f/ff/Data20140803.tar.gz")
set(DATA_MD5SUM "1f0016cf3b1cf31250293de08d7fe759")
message("Downloading data files from crpropa.desy.de ~ 30 MB")
if(CMAKE_VERSION GREATER 2.6)
	file(DOWNLOAD ${DATA_URL} ${CMAKE_SOURCE_DIR}/data.tar.gz STATUS status EXPECTED_MD5 ${DATA_MD5SUM})
else()
	execute_process(COMMAND wget ${DATA_URL} -O ${CMAKE_SOURCE_DIR}/data.tar.gz --no-check-certificate RESULT_VARIABLE status)
endif(CMAKE_VERSION GREATER 2.6)
message("Download complete: ${status}")
message("If the download failed, please manually obtain the file from " ${DATA_URL} ", rename it to data.tar.gz and place it in the source folder'")

file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/data)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/data.tar.gz WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data)

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
	src/module/OutputTXT.cpp
	src/module/OutputShell.cpp
	src/module/OutputROOT.cpp
	src/module/OutputCRPropa2.cpp
	src/module/PhotonOutput1D.cpp
	src/module/PhotonDINT.cpp
	src/module/PhotonDINT1D.cpp
	src/module/PhotonEleCa.cpp
	src/module/Tools.cpp
	src/magneticField/MagneticField.cpp
	src/magneticField/MagneticFieldGrid.cpp
	src/magneticField/TurbulentMagneticField.cpp
	src/magneticField/JF12Field.cpp
	${CRPROPA_EXTRA_SOURCES}
)
target_link_libraries(crpropa ${CRPROPA_EXTRA_LIBRARIES})

add_executable(crpropa-xmlrun src/main.cpp)
target_link_libraries(crpropa-xmlrun crpropa)

# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
add_definitions(-DCRPROPA_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS crpropa-xmlrun RUNTIME DESTINATION bin)
install(TARGETS crpropa DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY data/ DESTINATION share/crpropa/ PATTERN ".git" EXCLUDE)

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

	add_executable(testPropagation test/testPropagation.cpp)
	target_link_libraries(testPropagation crpropa gtest gtest_main pthread)
	add_test(testPropagation testPropagation)

	add_executable(testBreakCondition test/testBreakCondition.cpp)
	target_link_libraries(testBreakCondition crpropa gtest gtest_main pthread)
	add_test(testBreakCondition testBreakCondition)

	add_executable(testOutput test/testOutput.cpp)
	target_link_libraries(testOutput crpropa gtest gtest_main pthread)
	add_test(testOutput testOutput)

	add_executable(testInteraction test/testInteraction.cpp)
	target_link_libraries(testInteraction crpropa gtest gtest_main pthread)
	add_test(testInteraction testInteraction)

	add_executable(testSource test/testSource.cpp)
	target_link_libraries(testSource crpropa gtest gtest_main pthread)
	add_test(testSource testSource)

	if(WITH_GALACTIC_LENSES)
		include_directories(
			libs/GalacticMagneticLens/include
			)
		add_executable(testGalacticMagneticLens test/testMagneticLens.cpp)
		target_link_libraries(testGalacticMagneticLens crpropa gtest gtest_main pthread)
		add_test(testGalacticMagneticLens testGalacticMagneticLens)
	endif(WITH_GALACTIC_LENSES)

	if(ENABLE_PYTHON)
		CONFIGURE_FILE(test/testMagneticLensPythonInterface.py.in testMagneticLensPythonInterface.py) 
		ADD_TEST(testMagneticLensPythonInterface ${PYTHON_EXECUTABLE} testMagneticLensPythonInterface.py)
	endif(ENABLE_PYTHON)
endif(ENABLE_TESTING)

# ----------------------------------------------------------------------------
# Python
# ----------------------------------------------------------------------------
option(ENABLE_PYTHON "Create python library via SWIG" ON)
if(ENABLE_PYTHON)
	FIND_PACKAGE(SWIG REQUIRED)
	include(python/Python.cmake)
	include_directories(${PYTHON_INCLUDE_PATH})

	# tries to import numpy
	execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/checkNumpy.py" OUTPUT_VARIABLE numpyIncludePath)
	if(numpyIncludePath)
	  MESSAGE(STATUS "Found numpy headers in " ${numpyIncludePath})
		SET(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
		list(APPEND CRPROPA_SWIG_DEFINES -DWITHNUMPY)
	  include_directories(${numpyIncludePath})
	else()
	  MESSAGE(STATUS "Numpy not found.")
	endif()

	file(GLOB_RECURSE CRPROPA_INCLUDES include/*.h)
	set_source_files_properties( ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx PROPERTIES GENERATED true )
	add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx
		COMMAND swig -c++ -python -I${CMAKE_SOURCE_DIR}/include ${CRPROPA_SWIG_DEFINES} -o ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/python/crpropa.i
		DEPENDS ${CMAKE_SOURCE_DIR}/python/crpropa.i ${CRPROPA_INCLUDES} )

	add_library(crpropa-swig MODULE ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx)
	set_target_properties(crpropa-swig PROPERTIES PREFIX "")
	set_target_properties(crpropa-swig PROPERTIES OUTPUT_NAME "_crpropa")
	target_link_libraries(crpropa-swig crpropa ${PYTHON_LIBRARIES})

	install(DIRECTORY "${CMAKE_SOURCE_DIR}/python/crpropa" DESTINATION ${PYTHON_SITE_PACKAGES})
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/crpropa.py" DESTINATION ${PYTHON_SITE_PACKAGES}/crpropa)
	install(TARGETS crpropa-swig LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}/crpropa)
endif(ENABLE_PYTHON)
