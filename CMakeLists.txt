cmake_minimum_required(VERSION 2.6)
project(ModularPropagationCode)

find_program(GADGET_CONFIG_COMMAND "gadget-config") 
message("gadget-config: ${GADGET_CONFIG_COMMAND}")
execute_process(COMMAND ${GADGET_CONFIG_COMMAND} "--include" OUTPUT_VARIABLE GADGET_INCLUDES OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Gadget includes: ${GADGET_INCLUDES}")

execute_process(COMMAND ${GADGET_CONFIG_COMMAND} "--libs" OUTPUT_VARIABLE GADGET_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Gadget includes: ${GADGET_LIBS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GADGET_INCLUDES} ${OPTIONAL_CFLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GADGET_LIBS} ${OPTIONAL_LIBS}")

add_executable(mpc
    src/main.cpp
    src/ThreeVector.cpp
)

TARGET_LINK_LIBRARIES(mpc ${GADGET_LIBS})

# gtest
enable_testing()
find_package(GTest REQUIRED)
include_directories(src ${GTEST_INCLUDE_DIRS})

add_executable(testParticle test/testParticle.cpp)
target_link_libraries(testParticle ${GTEST_BOTH_LIBRARIES} pthread)
add_test(test1 testParticle)

add_executable(testERK test/testERK.cpp src/ThreeVector.cpp)
target_link_libraries(testERK ${GTEST_BOTH_LIBRARIES} pthread)
add_test(test2 testERK)