cmake_minimum_required(VERSION 3.14)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(HepPID STATIC
	src/ParticleIDMethods.cc
	src/ParticleName.cc
)

set_target_properties(HepPID PROPERTIES COMPILE_FLAGS -fPIC)
