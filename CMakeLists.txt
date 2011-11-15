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