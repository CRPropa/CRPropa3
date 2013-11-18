# SQLITE3_INCLUDE_DIR = path to Quimby directory
# SQLITE3_LIBRARY = libQuimby.a
# SQLITE3_FOUND = true if Quimby is found

find_path(SQLITE3_INCLUDE_DIR sqlite3.h)
find_library(SQLITE3_LIBRARY libsqlite3)

set(SQLITE3_FOUND FALSE)
if(SQLITE3_INCLUDE_DIR AND SQLITE3_LIBRARY)
    set(QUIMBY_FOUND TRUE)
    MESSAGE(STATUS "SQLite3: Found!")
else()
    MESSAGE(STATUS "SQLite3: NOT Found!")    
endif()

MESSAGE(STATUS "  Include:     ${SQLITE3_INCLUDE_DIR}")
MESSAGE(STATUS "  Library:     ${SQLITE3_LIBRARY}")

mark_as_advanced(SQLITE3_INCLUDE_DIR SQLITE3_LIBRARY QUIMBY_FOUND)
