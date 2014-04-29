# - Try to find muParser
# Once done this will define
#  MUPARSER_FOUND - System has LibXml2
#  MUPARSER_INCLUDE_DIRS - The LibXml2 include directories
#  MUPARSER_LIBRARIES - The libraries needed to use LibXml2
#  MUPARSER_DEFINITIONS - Compiler switches required for using LibXml2

find_package(PkgConfig)
pkg_check_modules(PC_MUPARSER QUIET muparser)
set(MUPARSER_DEFINITIONS ${PC_MUPARSER_CFLAGS_OTHER})

find_path(MUPARSER_INCLUDE_DIR muParser.h
          HINTS ${PC_MUPARSER_INCLUDEDIR} ${PC_MUPARSER_INCLUDE_DIRS})

find_library(MUPARSER_LIBRARY NAMES muparser libmuparser
             HINTS ${PC_MUPARSER_LIBDIR} ${PC_MUPARSER_LIBRARY_DIRS} )

set(MUPARSER_LIBRARIES ${MUPARSER_LIBRARY} )
set(MUPARSER_INCLUDE_DIRS ${MUPARSER_INCLUDE_DIR} )

set(MUPARSER_FOUND FALSE)
if(MUPARSER_INCLUDE_DIR AND MUPARSER_LIBRARY)
    set(MUPARSER_FOUND TRUE)
    MESSAGE(STATUS "muParser: Found!")
else()
    MESSAGE(STATUS "muParser: NOT Found!")
endif()

MESSAGE(STATUS "  Include:     ${MUPARSER_INCLUDE_DIR}")
MESSAGE(STATUS "  Library:     ${MUPARSER_LIBRARY}")

mark_as_advanced(MUPARSER_INCLUDE_DIR MUPARSER_LIBRARY )
