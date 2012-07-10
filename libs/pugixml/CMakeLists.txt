cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(pugixml
    pugixml
)

SET_TARGET_PROPERTIES( pugixml PROPERTIES COMPILE_FLAGS -fPIC)
