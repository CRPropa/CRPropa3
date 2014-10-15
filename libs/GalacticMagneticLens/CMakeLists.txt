cmake_minimum_required(VERSION 2.6)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${Boost_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/libs/HepPID/include
	)

add_library(galacticmagneticlens STATIC
	src/MagneticLens.cpp
	src/ModelMatrix.cpp
	src/Pixelization.cc
	src/Random.cc
	src/ParticleMapsContainer.cc
	)

target_link_libraries(galacticmagneticlens HepPID)
SET_TARGET_PROPERTIES(galacticmagneticlens PROPERTIES COMPILE_FLAGS -fPIC)

