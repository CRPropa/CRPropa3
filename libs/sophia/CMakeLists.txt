cmake_minimum_required(VERSION 2.6)

enable_language(Fortran)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(sophia STATIC
    sophia_interface.f
)

set_target_properties(sophia PROPERTIES COMPILE_FLAGS -fPIC)
