cmake_minimum_required(VERSION 2.6)
project(ModularPropagationCode)

include_directories(include)

# libmpc
add_library(mpc
    src/Candidate.cpp
    src/DeflectionCK.cpp
    src/Vector3.cpp
    src/XMLImport.cpp
    src/Module.cpp
    src/ParticleState.cpp
    src/ModuleChain.cpp
    src/SphMagneticField.cpp
    src/interaction/ElectronPairProduction.cpp
    src/interaction/Decay.cpp
    src/pugixml/pugixml.cpp
)
set(LIBS ${LIBS} mpc)

# components with external dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# fftw needed for TurbulentMagneticField
find_package(FFTW)
if (FFTW_FOUND)
    add_library(turbulentfield src/TurbulentMagneticField.cpp)
    set(LIBS ${LIBS} turbulentfield ${FFTW_LIBRARIES})
endif (FFTW_FOUND)

# freeglut and GLU needed for GlutDisplay
find_package(FreeGlut)
find_package(OpenGl)
if (Freeglut_FOUND)
    if (OPENGL_GLU_FOUND)
        add_library(glutdisplay src/GlutDisplay.cpp)
        set(LIBS ${LIBS} glutdisplay ${Freeglut_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
    endif (OPENGL_GLU_FOUND)
endif (Freeglut_FOUND)

# executable
add_executable(mpc-run src/main.cpp)
TARGET_LINK_LIBRARIES(mpc-run ${LIBS})

# testing
find_package(GTest)
if (GTEST_FOUND)
    enable_testing()

    add_executable(testCore test/testCore.cpp)
    target_link_libraries(testCore ${LIBS} ${GTEST_BOTH_LIBRARIES} pthread)
    add_test(test1 testCore)

    add_executable(testInteractions test/testInteractions.cpp)
    target_link_libraries(testInteractions ${LIBS} ${GTEST_BOTH_LIBRARIES} pthread)
    add_test(test2 testInteractions)
    
    add_executable(testMagneticField test/testMagneticField.cpp)
    target_link_libraries(testMagneticField ${LIBS} ${GTEST_BOTH_LIBRARIES} pthread)
    add_test(test3 testMagneticField)
endif (GTEST_FOUND)

#INCLUDE (Python.cmake)
#include_directories(${PYTHON_INCLUDE_PATH})

#FILE(GLOB_RECURSE MPC_INCLUDES include/*.h)
#message("Includes: ${MPC_INCLUDES}")
#SET_SOURCE_FILES_PROPERTIES( ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py PROPERTIES GENERATED true )
#ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py
#                COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/python/genbindings.py ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py ${CMAKE_SOURCE_DIR}/include ${MPC_INCLUDES}
#                DEPENDS ${CMAKE_SOURCE_DIR}/python/genbindings.py ${MPC_INCLUDES} )


#SET_SOURCE_FILES_PROPERTIES( ${CMAKE_CURRENT_BINARY_DIR}/mpc.cpp PROPERTIES GENERATED true )
#ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mpc.cpp
#                COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py > ${CMAKE_CURRENT_BINARY_DIR}/mpc.cpp
#                DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/genmodule.py )
                
#add_library(mpc_python MODULE ${CMAKE_CURRENT_BINARY_DIR}/mpc.cpp)
#set_target_properties(mpc_python PROPERTIES PREFIX "")
#set_target_properties(mpc_python PROPERTIES OUTPUT_NAME "mpc") 
#target_link_libraries(mpc_python libppropa ${PYTHON_LIBRARIES})