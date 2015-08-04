# SQLITE3_INCLUDE_DIR = path to SAGA directory
# SQLITE3_LIBRARY = libsaga.so
# SQLITE3_FOUND = true if SAGA is found

find_path(SQLITE3_INCLUDE_DIR sqlite3.h)
find_library(SQLITE3_LIBRARY libsqlite3)

set(SQLITE3_FOUND FALSE)
if(SQLITE3_INCLUDE_DIR AND SQLITE3_LIBRARY)
    set(SAGA_FOUND TRUE)
    MESSAGE(STATUS "SQLite3: Found!")
    include_directories(${SQLITE3_INCLUDE_DIR})
else()
    MESSAGE(STATUS "SQLite3: NOT Found!")    
endif()

MESSAGE(STATUS "  Include:     ${SQLITE3_INCLUDE_DIR}")
MESSAGE(STATUS "  Library:     ${SQLITE3_LIBRARY}")

mark_as_advanced(SQLITE3_INCLUDE_DIR SQLITE3_LIBRARY SAGA_FOUND)
