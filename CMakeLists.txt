cmake_minimum_required(VERSION 2.6)

project(ppropa CXX C)

include_directories(include)
#lib
add_library(libppropa
    src/ThreeVector.cpp
)

#exe
add_executable(ppropa
    src/main.cpp
)
target_link_libraries(ppropa libppropa)

# gtest
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# test 1
add_executable(testParticle test/testParticle.cpp)
target_link_libraries(testParticle libppropa ${GTEST_BOTH_LIBRARIES} pthread)
add_test(test1 testParticle)

#test 2
add_executable(testDeflectionCK test/testDeflectionCK.cpp)
target_link_libraries(testDeflectionCK libppropa ${GTEST_BOTH_LIBRARIES} pthread)
add_test(test2 testDeflectionCK)

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
