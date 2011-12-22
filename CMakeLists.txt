cmake_minimum_required(VERSION 2.6)
project(ModularPropagationCode)

include_directories(include)
# components with external dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set (MPC_EXTRA_SOURCES)
set (MPC_EXTRA_LIBRARIES)

# fftw needed for TurbulentMagneticField
find_package(FFTW)
if (FFTW_FOUND)
    list( APPEND MPC_EXTRA_SOURCES src/TurbulentMagneticField.cpp)
    list( APPEND MPC_EXTRA_LIBRARIES ${FFTW_LIBRARIES})
endif (FFTW_FOUND)

# freeglut and GLU needed for GlutDisplay
find_package(FreeGlut)
find_package(OpenGl)
if (Freeglut_FOUND)
    if (OPENGL_GLU_FOUND)
        list( APPEND MPC_EXTRA_SOURCES src/GlutDisplay.cpp)
        list( APPEND MPC_EXTRA_LIBRARIES ${Freeglut_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
    endif (OPENGL_GLU_FOUND)
endif (Freeglut_FOUND)

# libmpc
add_library(mpc SHARED
    src/Candidate.cpp
    src/DeflectionCK.cpp
    src/Vector3.cpp
    src/XMLImport.cpp
    src/Module.cpp
    src/ParticleState.cpp
    src/ModuleChain.cpp
    src/SphMagneticField.cpp
    src/interaction/ElectronPairProduction.cpp
    src/interaction/PhotoDisintegration.cpp
    src/interaction/Decay.cpp
    src/pugixml/pugixml.cpp
    ${MPC_EXTRA_SOURCES}
)
TARGET_LINK_LIBRARIES(mpc ${MPC_EXTRA_LIBRARIES})

# executable
add_executable(mpc-run src/main.cpp)
TARGET_LINK_LIBRARIES(mpc-run mpc)

# testing
find_package(GTest)
if (GTEST_FOUND)
    enable_testing()

    add_executable(testCore test/testCore.cpp)
    target_link_libraries(testCore mpc ${GTEST_BOTH_LIBRARIES} pthread)
    add_test(test1 testCore)

    add_executable(testInteractions test/testInteractions.cpp)
    target_link_libraries(testInteractions mpc ${GTEST_BOTH_LIBRARIES} pthread)
    add_test(test2 testInteractions)
    
    add_executable(testMagneticField test/testMagneticField.cpp)
    target_link_libraries(testMagneticField mpc ${GTEST_BOTH_LIBRARIES} pthread)
    add_test(test3 testMagneticField)
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