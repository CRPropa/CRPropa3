# SAGA_INCLUDE_DIR = path to SAGA directory
# SAGA_LIBRARY libsaga.a
# SAGA_FOUND = true if SAGA is found

find_path(SAGA_INCLUDE_DIR AMRgrid.h)
find_library(SAGA_LIBRARY libSAGA)

set(SAGA_FOUND FALSE)
if(SAGA_INCLUDE_DIR AND SAGA_LIBRARY)
    set(SAGA_FOUND TRUE)
    MESSAGE(STATUS "SAGA: Found!")
else()
    MESSAGE(STATUS "SAGA: NOT Found!")    
endif()

MESSAGE(STATUS "  Include:     ${SAGA_INCLUDE_DIR}")
MESSAGE(STATUS "  Library:     ${SAGA_LIBRARY}")

mark_as_advanced(SAGA_INCLUDE_DIR SAGA_LIBRARY SAGA_FOUND)
