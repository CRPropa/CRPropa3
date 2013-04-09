cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(dint STATIC
    src/advance.cpp
    src/background.cpp
    src/binfread.cpp
    src/check.cpp
    src/cvector.cpp
    src/decay.cpp
    src/deriv.cpp
    src/error.cpp
    src/final.cpp
    src/fold.cpp
    src/frag.cpp
    src/gauleg.cpp
    src/inject.cpp
    src/io_util.cpp
    src/load.cpp
    src/math_util.cpp
    src/prepare.cpp
    src/prop_second.cpp
    src/rate.cpp
    src/spectrum.cpp
    src/sync.cpp
    src/vector.cpp
)

SET_TARGET_PROPERTIES(dint PROPERTIES COMPILE_FLAGS -fPIC)

