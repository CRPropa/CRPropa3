cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(HepPID STATIC
	src/ParticleIDMethods.cc
	src/ParticleName.cc
	src/Version.cc
	src/translateEvtGen.cc
	src/translateGeanttoPDT.cc
	src/translateHerwig.cc
	src/translateIsajet.cc
	src/translatePDG.cc
	src/translatePDTtoGeant.cc
	src/translatePythia.cc
	src/translateQQ.cc
)

SET_TARGET_PROPERTIES(HepPID PROPERTIES COMPILE_FLAGS -fPIC)
