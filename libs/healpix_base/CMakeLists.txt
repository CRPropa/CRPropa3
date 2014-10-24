cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(healpix_base STATIC
	error_handling.cc
	geom_utils.cc
	healpix_base.cc
	healpix_tables.cc
	pointing.cc
	)

SET_TARGET_PROPERTIES( healpix_base PROPERTIES COMPILE_FLAGS -fPIC)
