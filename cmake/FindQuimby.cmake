# QUIMBY_INCLUDE_DIR = path to Quimby directory
# QUIMBY_LIBRARY = libQuimby.a
# QUIMBY_FOUND = true if Quimby is found

find_path(QUIMBY_INCLUDE_DIR quimby/SmoothParticle.h)
find_library(QUIMBY_LIBRARY Quimby)

set(QUIMBY_FOUND FALSE)
if(QUIMBY_INCLUDE_DIR AND QUIMBY_LIBRARY)
    set(QUIMBY_FOUND TRUE)
    MESSAGE(STATUS "Quimby: Found!")
else()
    MESSAGE(STATUS "Quimby: NOT Found!")    
endif()

MESSAGE(STATUS "  Include:     ${QUIMBY_INCLUDE_DIR}")
MESSAGE(STATUS "  Library:     ${QUIMBY_LIBRARY}")

mark_as_advanced(QUIMBY_INCLUDE_DIR QUIMBY_LIBRARY QUIMBY_FOUND)
