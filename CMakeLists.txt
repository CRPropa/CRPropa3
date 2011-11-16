cmake_minimum_required(VERSION 2.6)

project(mpc CXX C)

add_executable(mpc
    src/main.cpp
    src/ThreeVector.cpp
)

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