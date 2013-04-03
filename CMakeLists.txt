cmake_minimum_required(VERSION 2.6)
project(ModularPropagationCode Fortran C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(MPC_EXTRA_SOURCES)
set(MPC_EXTRA_INCLUDES)
set(MPC_EXTRA_LIBRARIES)
set(MPC_SWIG_DEFINES)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
# googletest (provided)
option(ENABLE_TESTING "Build tests and enable test target" ON)
if(ENABLE_TESTING)
	include_directories(libs/gtest/include)
	add_subdirectory(libs/gtest)
endif(ENABLE_TESTING)

# kiss (provided)
add_subdirectory(libs/kiss)
list(APPEND MPC_EXTRA_LIBRARIES kiss)
list(APPEND MPC_EXTRA_INCLUDES libs/kiss/include)

# pugixml (provided)
add_subdirectory(libs/pugixml)
list(APPEND MPC_EXTRA_LIBRARIES pugixml)
list(APPEND MPC_EXTRA_INCLUDES libs/pugixml)

# HepID (provided)
add_subdirectory(libs/HepPID)
list(APPEND MPC_EXTRA_LIBRARIES HepPID)
list(APPEND MPC_EXTRA_INCLUDES libs/HepPID/include)

# SOPHIA (provided)
add_subdirectory(libs/sophia)
list(APPEND MPC_EXTRA_LIBRARIES sophia gfortran)
list(APPEND MPC_EXTRA_INCLUDES libs/sophia)

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
        list(APPEND MPC_EXTRA_INCLUDES ${FFTW3F_INCLUDE_DIR})
        list(APPEND MPC_EXTRA_LIBRARIES ${FFTW3F_LIBRARY})
        add_definitions(-DMPC_HAVE_FFTW3F)
        list(APPEND MPC_SWIG_DEFINES -DMPC_HAVE_FFTW3F)
    endif(FFTW3F_FOUND)
endif(ENABLE_FFTW3F)

# gadget (optional for SPH magnetic fields)
find_package(Gadget)
if(GADGET_FOUND)
    list(APPEND MPC_EXTRA_INCLUDES ${GADGET_INCLUDE_DIR})
    list(APPEND MPC_EXTRA_LIBRARIES ${GADGET_LIBRARY})
    add_definitions (-DMPC_HAVE_GADGET)
    list(APPEND MPC_SWIG_DEFINES -DMPC_HAVE_GADGET)
    list(APPEND MPC_SWIG_DEFINES -I${GADGET_INCLUDE_DIR}/../share/gadget)
    list(APPEND MPC_SWIG_DEFINES -I${GADGET_INCLUDE_DIR})
endif(GADGET_FOUND)

# ROOT (optional for ROOT output)
option(ENABLE_ROOT "ROOT Output" ON)
if(ENABLE_ROOT)
    find_package(ROOT)
    if(ROOT_FOUND)
    	FIND_PROGRAM(ROOT_CONFIG root-config REQUIRED ROOT_CONFIG_PATH) 
    	execute_process(COMMAND ${ROOT_CONFIG} "--cflags" OUTPUT_VARIABLE ROOT_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
    	execute_process(COMMAND ${ROOT_CONFIG} "--libs" OUTPUT_VARIABLE ROOT_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
    	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CFLAGS}")
    	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ROOT_LIBS}")
    	set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ROOT_LIBS}")
    	add_definitions(-DMPC_HAVE_ROOT)
    	list(APPEND MPC_SWIG_DEFINES -DMPC_HAVE_ROOT)
    endif(ROOT_FOUND)
endif(ENABLE_ROOT)

# Google Performance Tools (optional as possible performance tweak)
find_package(GooglePerfTools)
set(TCMALLOC)
if(GOOGLE_PERFTOOLS_FOUND AND NOT NO_TCMALLOC)
    set(TCMALLOC ${TCMALLOC_LIBRARY})
endif(GOOGLE_PERFTOOLS_FOUND AND NOT NO_TCMALLOC)

# ----------------------------------------------------------------------------
# Library and Binary
# ----------------------------------------------------------------------------
include_directories(include ${MPC_EXTRA_INCLUDES})
add_library(mpc SHARED
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
    src/module/StochasticInteraction.cpp
    src/module/NuclearDecay.cpp
    src/module/PhotoPionProduction.cpp
    src/module/PhotoDisintegration.cpp
    src/module/Redshift.cpp
    src/module/Output.cpp
    src/module/OutputROOT.cpp
    src/module/OutputCRPropa2.cpp
    src/module/Tools.cpp
    src/magneticField/MagneticField.cpp
    src/magneticField/MagneticFieldGrid.cpp
    src/magneticField/TurbulentMagneticField.cpp
    src/magneticField/JF2012Field.cpp
    ${MPC_EXTRA_SOURCES}
)
target_link_libraries(mpc ${MPC_EXTRA_LIBRARIES})

add_executable(mpc-run src/main.cpp)
target_link_libraries(mpc-run mpc)

# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
add_definitions(-DMPC_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS mpc-run RUNTIME DESTINATION bin)
install(TARGETS mpc DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY data/ DESTINATION share/mpc/)

# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------
if(ENABLE_TESTING)
    enable_testing()
	add_executable(testCore test/testCore.cpp)
	target_link_libraries(testCore mpc gtest gtest_main pthread)
	add_test(testCore testCore)
	
	add_executable(testVector3 test/testVector3.cpp)
	target_link_libraries(testVector3 mpc gtest gtest_main pthread)
	add_test(testVector3 testVector3)
	
	add_executable(testModuleList test/testModuleList.cpp)
	target_link_libraries(testModuleList mpc gtest gtest_main pthread)
	add_test(testModuleList testModuleList)

	add_executable(testMagneticField test/testMagneticField.cpp)
	target_link_libraries(testMagneticField mpc gtest gtest_main pthread)
	add_test(testMagneticField testMagneticField)

	add_executable(testPropagation test/testPropagation.cpp)
	target_link_libraries(testPropagation mpc gtest gtest_main pthread)
	add_test(testPropagation testPropagation)
     
	add_executable(testBreakCondition test/testBreakCondition.cpp)
	target_link_libraries(testBreakCondition mpc gtest gtest_main pthread)
	add_test(testBreakCondition testBreakCondition)
    
    add_executable(testOutput test/testOutput.cpp)
	target_link_libraries(testOutput mpc gtest gtest_main pthread)
	add_test(testOutput testOutput)
    
	add_executable(testInteraction test/testInteraction.cpp)
	target_link_libraries(testInteraction mpc gtest gtest_main pthread)
	add_test(testInteraction testInteraction)
	
	add_executable(testSource test/testSource.cpp)
	target_link_libraries(testSource mpc gtest gtest_main pthread)
	add_test(testSource testSource)
endif(ENABLE_TESTING)

# ----------------------------------------------------------------------------
# Python
# ----------------------------------------------------------------------------
option(ENABLE_PYTHON "Create python library via SWIG" ON)
if(ENABLE_PYTHON)
	include(python/Python.cmake)
	include_directories(${PYTHON_INCLUDE_PATH})
	
	file(GLOB_RECURSE MPC_INCLUDES include/*.h)
	set_source_files_properties( ${CMAKE_CURRENT_BINARY_DIR}/mpc_wrap.cxx PROPERTIES GENERATED true )
	add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mpc_wrap.cxx
	                COMMAND swig -c++ -python ${MPC_SWIG_DEFINES} -I${CMAKE_SOURCE_DIR}/include -o ${CMAKE_CURRENT_BINARY_DIR}/mpc_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/python/mpc.i 
	                DEPENDS ${CMAKE_SOURCE_DIR}/python/mpc.i  ${MPC_INCLUDES} )
	add_library(mpc-swig MODULE ${CMAKE_CURRENT_BINARY_DIR}/mpc_wrap.cxx)
	set_target_properties(mpc-swig PROPERTIES PREFIX "")
	set_target_properties(mpc-swig PROPERTIES OUTPUT_NAME "_mpc") 
	target_link_libraries(mpc-swig mpc ${PYTHON_LIBRARIES})
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/mpc.py" DESTINATION ${PYTHON_SITE_PACKAGES})
	install(TARGETS mpc-swig LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES})
endif(ENABLE_PYTHON)
