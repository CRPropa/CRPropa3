cmake_minimum_required(VERSION 2.6)
project(ModularPropagationCode)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}" ${CMAKE_MODULE_PATH})

set (MPC_EXTRA_SOURCES)
set (MPC_EXTRA_INCLUDES)
set (MPC_EXTRA_LIBRARIES)

# GSL required
find_package(GSL REQUIRED)
list (APPEND MPC_EXTRA_LIBRARIES ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

# fftw needed for TurbulentMagneticField
find_package(FFTW)
if (FFTW_FOUND)
    list( APPEND MPC_EXTRA_SOURCES src/magneticField/turbulentMagneticFieldGrid.cpp)
    list( APPEND MPC_EXTRA_LIBRARIES ${FFTW_LIBRARY})
endif (FFTW_FOUND)

# gadget needed for SphMagneticField
find_package(Gadget)
if (GADGET_FOUND)
    list( APPEND MPC_EXTRA_SOURCES src/magneticField/sphMagneticField.cpp)
    list( APPEND MPC_EXTRA_INCLUDES ${GADGET_INCLUDE_DIR})
    list( APPEND MPC_EXTRA_LIBRARIES ${GADGET_LIBRARY})
endif (GADGET_FOUND)

# freeglut and GLU needed for GlutDisplay
find_package(FreeGlut)
find_package(OpenGl)
if (Freeglut_FOUND)
    if (OPENGL_GLU_FOUND)
        list( APPEND MPC_EXTRA_SOURCES src/module/GlutDisplay.cpp)
        list( APPEND MPC_EXTRA_LIBRARIES ${Freeglut_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
    endif (OPENGL_GLU_FOUND)
endif (Freeglut_FOUND)

# pugixml provided
add_library (pugixml libs/pugixml/pugixml.cpp)
SET_TARGET_PROPERTIES( pugixml PROPERTIES COMPILE_FLAGS -fPIC)
list (APPEND MPC_EXTRA_LIBRARIES pugixml)
list (APPEND MPC_EXTRA_INCLUDES libs/pugixml)

# kiss provided
add_subdirectory(libs/kiss)
list (APPEND MPC_EXTRA_LIBRARIES kiss)
list (APPEND MPC_EXTRA_INCLUDES libs/kiss/include)

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(include ${MPC_EXTRA_INCLUDES})

add_library(mpc SHARED
    src/Random.cpp
    src/Vector3.cpp
    src/XMLImport.cpp
    src/ModuleChain.cpp
    src/ModuleList.cpp
    src/Module.cpp
    src/Candidate.cpp
    src/IO.cpp
    src/ParticleState.cpp
    src/module/common.cpp
    src/module/SimplePropagation.cpp
    src/module/DeflectionCK.cpp
    src/module/ElectronPairProduction.cpp
    src/module/PhotoPionProduction.cpp
    src/module/PhotoDisintegration.cpp
    src/module/NuclearDecay.cpp
    src/module/Redshift.cpp
    src/magneticField/magneticFieldGrid.cpp
    ${MPC_EXTRA_SOURCES}
)
TARGET_LINK_LIBRARIES(mpc ${MPC_EXTRA_LIBRARIES})

add_executable(mpc-run src/main.cpp)
TARGET_LINK_LIBRARIES(mpc-run mpc)

add_executable(openmp example/openmp.cpp)
TARGET_LINK_LIBRARIES(openmp mpc)

# MPI needed for mpc-mpi
find_package(MPI)
if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})

    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

    add_executable(mpc-mpi mpi/main mpi/Slave mpi/Master mpi/Job)
    TARGET_LINK_LIBRARIES(mpc-mpi mpc ${MPI_LIBRARIES})
endif()


# GTest needed for testing
find_package(GTest)
if (GTEST_FOUND)
    enable_testing()
    
    add_executable(testCommon test/testCommon.cpp)
    target_link_libraries(testCommon mpc ${GTEST_BOTH_LIBRARIES} pthread)

    add_executable(testCore test/testCore.cpp)
    target_link_libraries(testCore mpc ${GTEST_BOTH_LIBRARIES} pthread)

    add_executable(testInteractions test/testInteractions.cpp)
    target_link_libraries(testInteractions mpc ${GTEST_BOTH_LIBRARIES} pthread)
    
    add_executable(testMagneticField test/testMagneticField.cpp)
    target_link_libraries(testMagneticField mpc ${GTEST_BOTH_LIBRARIES} pthread)
    
    add_executable(testBreakCondition test/testBreakCondition.cpp)
    target_link_libraries(testBreakCondition mpc ${GTEST_BOTH_LIBRARIES} pthread)
endif (GTEST_FOUND)

INCLUDE (Python.cmake)
include_directories(${PYTHON_INCLUDE_PATH})

FILE(GLOB_RECURSE MPC_INCLUDES include/*.h)
#message("Includes: ${MPC_INCLUDES}")
#SET_SOURCE_FILES_PROPERTIES( ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py PROPERTIES GENERATED true )
#ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py
#                COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/python/genbindings.py ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py ${CMAKE_SOURCE_DIR}/include ${MPC_INCLUDES}
#                DEPENDS ${CMAKE_SOURCE_DIR}/python/genbindings.py ${MPC_INCLUDES} )


#SET_SOURCE_FILES_PROPERTIES( ${CMAKE_CURRENT_BINARY_DIR}/mpc_python.cpp PROPERTIES GENERATED true )
#ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mpc_python.cpp
#                COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py > ${CMAKE_CURRENT_BINARY_DIR}/mpc_python.cpp
#                DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py )
                
#add_library(mpc_python MODULE ${CMAKE_CURRENT_BINARY_DIR}/mpc_python.cpp)
#set_target_properties(mpc_python PROPERTIES PREFIX "")
#set_target_properties(mpc_python PROPERTIES OUTPUT_NAME "mpc") 
#target_link_libraries(mpc_python libppropa ${PYTHON_LIBRARIES})

SET_SOURCE_FILES_PROPERTIES( ${CMAKE_CURRENT_BINARY_DIR}/mpc_wrap.cxx PROPERTIES GENERATED true )
ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mpc_wrap.cxx
                COMMAND swig -c++ -python -I${CMAKE_SOURCE_DIR}/include -o ${CMAKE_CURRENT_BINARY_DIR}/mpc_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/python/mpc.i 
                DEPENDS ${CMAKE_SOURCE_DIR}/python/mpc.i  ${MPC_INCLUDES} )
add_library(mpc-swig MODULE ${CMAKE_CURRENT_BINARY_DIR}/mpc_wrap.cxx)
set_target_properties(mpc-swig PROPERTIES PREFIX "")
set_target_properties(mpc-swig PROPERTIES OUTPUT_NAME "_mpc") 
target_link_libraries(mpc-swig mpc ${PYTHON_LIBRARIES})