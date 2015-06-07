cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(eleca STATIC
    src/Common
    src/EnergyLoss
    src/Particle
    src/Process
    src/Propagation
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET_TARGET_PROPERTIES(eleca PROPERTIES COMPILE_FLAGS "-fPIC -ffast-math -D__extern_always_inline=inline")
else()
    SET_TARGET_PROPERTIES(eleca PROPERTIES COMPILE_FLAGS "-fPIC -ffast-math")
endif()
