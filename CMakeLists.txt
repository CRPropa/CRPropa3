cmake_minimum_required(VERSION 2.6)
project(ModularPropagationCode)

include_directories(include)
#lib
add_library(libmpc
    src/ThreeVector.cpp
)

find_program(GADGET_CONFIG_COMMAND "gadget-config") 
message("gadget-config: ${GADGET_CONFIG_COMMAND}")
execute_process(COMMAND ${GADGET_CONFIG_COMMAND} "--include" OUTPUT_VARIABLE GADGET_INCLUDES OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Gadget includes: ${GADGET_INCLUDES}")

execute_process(COMMAND ${GADGET_CONFIG_COMMAND} "--libs" OUTPUT_VARIABLE GADGET_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Gadget includes: ${GADGET_LIBS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GADGET_INCLUDES} ${OPTIONAL_CFLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GADGET_LIBS} ${OPTIONAL_LIBS}")

# exe
add_executable(mpc
    src/main.cpp
)

TARGET_LINK_LIBRARIES(mpc ${GADGET_LIBS} libmpc glut GLU)

# gtest
enable_testing()
find_package(GTest REQUIRED)
include_directories(src ${GTEST_INCLUDE_DIRS})

# test 1
add_executable(testParticle test/testParticle.cpp)
target_link_libraries(testParticle ${GTEST_BOTH_LIBRARIES} pthread)
add_test(test1 testParticle)

# test 2
add_executable(testERK test/testERK.cpp src/ThreeVector.cpp)
target_link_libraries(testERK ${GTEST_BOTH_LIBRARIES} pthread)
add_test(test2 testERK)


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