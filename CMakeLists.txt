cmake_minimum_required(VERSION 2.6)
project(CRPropa Fortran C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CRPROPA_EXTRA_SOURCES)
set(CRPROPA_EXTRA_INCLUDES)
set(CRPROPA_EXTRA_LIBRARIES)
set(CRPROPA_SWIG_DEFINES)

if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
endif(CMAKE_COMPILER_IS_GNUCXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
# googletest (provided)
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
list(APPEND CRPROPA_EXTRA_INCLUDES libs/EleCa)

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

# fftw3f (optional for turbulent magnetic fields)
option(ENABLE_FFTW3F "FFTW3F to create turbulent fields" ON)
if(ENABLE_FFTW3F)
    find_package(FFTW3F)
    if(FFTW3F_FOUND)
        list(APPEND CRPROPA_EXTRA_INCLUDES ${FFTW3F_INCLUDE_DIR})
        list(APPEND CRPROPA_EXTRA_LIBRARIES ${FFTW3F_LIBRARY})
        add_definitions(-DCRPROPA_HAVE_FFTW3F)
        list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_FFTW3F)
    endif(FFTW3F_FOUND)
endif(ENABLE_FFTW3F)

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

find_package(muParser)
if(MUPARSER_FOUND)
    list(APPEND CRPROPA_EXTRA_INCLUDES ${MUPARSER_INCLUDE_DIR})
    list(APPEND CRPROPA_EXTRA_LIBRARIES ${MUPARSER_LIBRARY})
    add_definitions (-DCRPROPA_HAVE_MUPARSER)
    list(APPEND CRPROPA_SWIG_DEFINES -DCRPROPA_HAVE_MUPARSER)
    list(APPEND CRPROPA_SWIG_DEFINES -I${MUPARSER_INCLUDE_DIR})
endif(MUPARSER_FOUND)


# Google Performance Tools (optional as possible performance tweak)
find_package(GooglePerfTools)
set(TCMALLOC)
if(GOOGLE_PERFTOOLS_FOUND AND NOT NO_TCMALLOC)
    set(TCMALLOC ${TCMALLOC_LIBRARY})
    list(APPEND CRPROPA_EXTRA_LIBRARIES profiler)
endif(GOOGLE_PERFTOOLS_FOUND AND NOT NO_TCMALLOC)

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
    src/Cosmology.cpp
    src/Source.cpp
    src/Common.cpp
    src/PhotonBackground.cpp
    src/GridTools.cpp
    src/XmlExecute.cpp
    src/module/BreakCondition.cpp
    src/module/Boundary.cpp
    src/module/Observer.cpp
    src/module/SimplePropagation.cpp
    src/module/DeflectionCK.cpp
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
endif(ENABLE_TESTING)

# ----------------------------------------------------------------------------
# Python
# ----------------------------------------------------------------------------
option(ENABLE_PYTHON "Create python library via SWIG" ON)
if(ENABLE_PYTHON)
	include(python/Python.cmake)
	include_directories(${PYTHON_INCLUDE_PATH})
	
	file(GLOB_RECURSE CRPROPA_INCLUDES include/*.h)
	set_source_files_properties( ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx PROPERTIES GENERATED true )
	add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx
	                COMMAND swig -c++ -python -I${CMAKE_SOURCE_DIR}/include ${CRPROPA_SWIG_DEFINES} -o ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/python/crpropa.i 
	                DEPENDS ${CMAKE_SOURCE_DIR}/python/crpropa.i ${CRPROPA_INCLUDES} )
	
	add_library(crpropa-swig MODULE ${CMAKE_CURRENT_BINARY_DIR}/crpropa_wrap.cxx)
	set_target_properties(crpropa-swig PROPERTIES PREFIX "")
	set_target_properties(crpropa-swig PROPERTIES OUTPUT_NAME "_crpropa") 
	target_link_libraries(crpropa-swig crpropa ${PYTHON_LIBRARIES})
	
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/crpropa.py" DESTINATION ${PYTHON_SITE_PACKAGES})
	install(TARGETS crpropa-swig LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES})
endif(ENABLE_PYTHON)
