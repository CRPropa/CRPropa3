cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(eleca STATIC
    src/Common
    src/EnergyLoss
    src/Particle
    src/Process
    src/Propagation
)

SET_TARGET_PROPERTIES(eleca PROPERTIES COMPILE_FLAGS "-fPIC -ffast-math")
