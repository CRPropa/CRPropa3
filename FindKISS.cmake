# KISS_INCLUDE_DIR = path to include directory
# KISS_LIBRARY = path to library
# KISS_FOUND = true if package found is found

find_path(KISS_INCLUDE_DIR kiss/convert.h)
find_library(KISS_LIBRARY kiss)

set(KISS_FOUND NOTFOUND)
if(KISS_INCLUDE_DIR AND KISS_LIBRARY)
    set(KISS_FOUND TRUE)
    MESSAGE(STATUS "Found kiss")
else()
    MESSAGE(STATUS "kiss NOT found")
endif()

MESSAGE(STATUS "  Include: ${KISS_INCLUDE_DIR}")
MESSAGE(STATUS "  Library: ${KISS_LIBRARY}")

mark_as_advanced(KISS_INCLUDE_DIR KISS_LIBRARY)
