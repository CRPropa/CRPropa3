###############################################################################
# FindCRPropa.cmake
# -----------------------------------------------------------------------------
# Locate an installed CRPropa (headers, library, SWIG interface)
#
# This is a standard CMake find-module: it is used through
#   find_package(CRPropa [REQUIRED] [QUIET] [<version>])
# because the file lives in a directory that is on CMAKE_MODULE_PATH.
#
# Every value can be controlled manually.
# If you set one of the output variables below (e.g. with -DCRPropa_INCLUDE_DIR=/path or in the cache), the module keeps your value and only auto-detects the ones left empty.
#
# Input hints (all optional)
# - `CRPropa_INSTALL_PREFIX`: CRPropa install prefix (cache/-D variable)
# - `CRPropa_DIR` (environment) alternative install-prefix hint
# - `CRPropa_ROOT` (cache or environment) is also set automatically through the find_* commands under policy CMP0074 (CMake >= 3.12).
# - `Python_EXECUTABLE`: used as a last-resort locator (set by find_package(Python))
#
# Output variables
# - `CRPropa_FOUND`: whether all required pieces were located
# - `CRPropa_INCLUDE_DIR`: directory containing CRPropa.h
# - `CRPropa_LIBRARY`: the crpropa shared library
# - `CRPropa_LIBRARY_DIR`: directory containing the crpropa shared library
# - `CRPropa_SWIG_PATH`: directory containing the SWIG interface files
# - `CRPropa_SWIG_INTERFACE_FILE`: full path to the selected SWIG interface file
# - `CRPropa_VERSION`: detected version (best effort; may be empty)
#
# Imported target (preferred way to call CRPropa from plugins):
#   CRPropa::CRPropa  
# This carries the library and its include directory. 
# Therefore, `target_link_libraries(foo CRPropa::CRPropa)` is all that is needed.
#
###############################################################################



# Expose the install prefix as a cache entry so it can be passed with -D and edited in the cache/GUI.
# Existing (user-provided) values are preserved.
set(CRPropa_INSTALL_PREFIX "${CRPropa_INSTALL_PREFIX}" CACHE PATH "CRPropa install prefix (leave empty to auto-detect)")

# The SWIG interface file to use depends on the builtin mode.
# Builtin is the default because, for unknown reasons, its absence causes segmentation faults on various systems.
option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" On)
if(ENABLE_SWIG_BUILTIN)
	set(CRPropa_SWIG_FILE "crpropa-builtin.i")
else()
	set(CRPropa_SWIG_FILE "crpropa.i")
endif(ENABLE_SWIG_BUILTIN)


set(CRPropa_HINTS
	${CRPropa_INSTALL_PREFIX}
	${CRPropa_INSTALL_PREFIX}/build
	$ENV{CRPropa_DIR}
	$ENV{CRPropa_DIR}/build
	# ${CMAKE_CURRENT_SOURCE_DIR}/../build
	# ${CMAKE_CURRENT_SOURCE_DIR}/..
)

# Common system-wide installation prefixes.
set(CRPropa_SYSTEM_PATHS
	/usr
	/usr/local
	/opt
	/opt/local
	/usr/local/opt
)
if(APPLE)
	list(APPEND CRPropa_SYSTEM_PATHS "$ENV{HOME}/.local")
endif()


# Note that `find_path` and `find_library` do not work when their cache variable is already set.
# Therefore, a manually provided CRPropa_INCLUDE_DIR / CRPropa_LIBRARY always wins over the auto-detection below.
find_path(CRPropa_INCLUDE_DIR
	NAMES CRPropa.h
	HINTS ${CRPropa_HINTS}
	PATHS ${CRPropa_SYSTEM_PATHS}
	PATH_SUFFIXES include
)
find_library(CRPropa_LIBRARY
	NAMES crpropa
	HINTS ${CRPropa_HINTS}
	PATHS ${CRPropa_SYSTEM_PATHS}
	PATH_SUFFIXES lib lib64
)
if(CRPropa_LIBRARY)
	get_filename_component(CRPropa_LIBRARY_DIR "${CRPropa_LIBRARY}" DIRECTORY)
endif(CRPropa_LIBRARY)

# SWIG interface directory
find_path(CRPropa_SWIG_PATH
	NAMES crpropa.i
	HINTS ${CRPropa_HINTS}
	PATHS ${CRPropa_SYSTEM_PATHS}
	PATH_SUFFIXES share/crpropa/swig_interface
)

# Locate the Python helper relative to this module file.
set(_findCRPropa "${CMAKE_CURRENT_LIST_DIR}/python/findCRPropa.py")


# If anything is still missing, ask the installed  CRPropa Python module where it lives (see python/findCRPropa.py).
# This is useful if there are multiple Python installations and the one used to build CRPropa is not the one CMake finds by default.
if(Python_EXECUTABLE AND EXISTS "${_findCRPropa}" AND (NOT CRPropa_SWIG_PATH OR NOT CRPropa_INCLUDE_DIR OR NOT CRPropa_LIBRARY))
	if(NOT CRPropa_SWIG_PATH)
		execute_process(
			COMMAND ${Python_EXECUTABLE} ${_findCRPropa} swig_interface
			OUTPUT_VARIABLE _crpropaSwig
			OUTPUT_STRIP_TRAILING_WHITESPACE
			RESULT_VARIABLE rcCRPropa
		)
		if(rcCRPropa EQUAL 0 AND _crpropaSwig)
			set(CRPropa_SWIG_PATH "${_crpropaSwig}" CACHE PATH "CRPropa SWIG interface directory" FORCE)
		endif()
	endif(NOT CRPropa_SWIG_PATH)

	if(NOT CRPropa_INCLUDE_DIR OR NOT CRPropa_LIBRARY)
		execute_process(
			COMMAND ${Python_EXECUTABLE} ${_findCRPropa} install_prefix
			OUTPUT_VARIABLE prefixCRPropa
			OUTPUT_STRIP_TRAILING_WHITESPACE
			RESULT_VARIABLE rcCRPropa
		)
		if(rcCRPropa EQUAL 0 AND prefixCRPropa)
			find_path(CRPropa_INCLUDE_DIR NAMES CRPropa.h
				HINTS ${prefixCRPropa}/include)
			find_library(CRPropa_LIBRARY NAMES crpropa
				HINTS ${prefixCRPropa}/lib ${prefixCRPropa}/lib64)
		endif()
	endif(NOT CRPropa_INCLUDE_DIR OR NOT CRPropa_LIBRARY)
endif()

# Try to detect version
if(Python_EXECUTABLE AND EXISTS "${_findCRPropa}" AND NOT CRPropa_VERSION)
	execute_process(
		COMMAND ${Python_EXECUTABLE} ${_findCRPropa} version
		OUTPUT_VARIABLE versionCRPropa
		OUTPUT_STRIP_TRAILING_WHITESPACE
		ERROR_QUIET
		RESULT_VARIABLE rcCRPropa
	)
	if(rcCRPropa EQUAL 0 AND versionCRPropa MATCHES "^[0-9]")
		set(CRPropa_VERSION "${versionCRPropa}")
	endif()
endif()

# Full path to the selected SWIG interface file, verified to exist.
# Note: the directory above is located via crpropa.i, but crpropa-builtin.i may also be used.
if(CRPropa_SWIG_PATH)
	set(CRPropa_SWIG_INTERFACE_FILE "${CRPropa_SWIG_PATH}/${CRPropa_SWIG_FILE}")
	if(NOT EXISTS "${CRPropa_SWIG_INTERFACE_FILE}")
		message(WARNING
			"Selected SWIG interface file was not found:\n"
			"  ${CRPropa_SWIG_INTERFACE_FILE}\n"
			"Toggle ENABLE_SWIG_BUILTIN or check your CRPropa installation.")
	endif(NOT EXISTS "${CRPropa_SWIG_INTERFACE_FILE}")
endif(CRPropa_SWIG_PATH)


# If something is still missing, give instructions (including manual flags).
if((NOT CRPropa_INCLUDE_DIR OR NOT CRPropa_LIBRARY OR NOT CRPropa_SWIG_PATH) AND NOT CRPropa_FIND_QUIETLY)
	message(STATUS "")
	message(STATUS "CRPropa could not be located automatically. You can set it manually with one of:")
	message(STATUS "  cmake -DCRPropa_INSTALL_PREFIX=/path/to/crpropa/install ..")
	message(STATUS "  cmake -DCRPropa_DIR=/path/to/crpropa/install ..  (or export CRPropa_DIR)")
	message(STATUS "or provide the individual paths:")
	message(STATUS "  -DCRPropa_INCLUDE_DIR=/path/to/include  (contains CRPropa.h)")
	message(STATUS "  -DCRPropa_LIBRARY=/path/to/lib/libcrpropa.so|.dylib")
	message(STATUS "  -DCRPropa_SWIG_PATH=/path/to/share/crpropa/swig_interface")
	message(STATUS "If CRPropa is installed as a Python module, make sure the interpreter")
	message(STATUS "  '${Python_EXECUTABLE}' can 'import crpropa' (or pass -DPython_EXECUTABLE=...).")
	message(STATUS "")
endif((NOT CRPropa_INCLUDE_DIR OR NOT CRPropa_LIBRARY OR NOT CRPropa_SWIG_PATH) AND NOT CRPropa_FIND_QUIETLY)

# Sets the variable CRPropa_FOUND.
include(FindPackageHandleStandardArgs)
set(versionArgsCRPropa)
if(CRPropa_VERSION)
	list(APPEND versionArgsCRPropa VERSION_VAR CRPropa_VERSION)
endif()
find_package_handle_standard_args(CRPropa
	REQUIRED_VARS
		CRPropa_INCLUDE_DIR
		CRPropa_LIBRARY
		CRPropa_SWIG_PATH
	${versionArgsCRPropa}
)

# Provide imported target for linking against; inherit the include directory transitively. 
# UNKNOWN IMPORTED is used because `find_library` does not tell us whether the library is static or shared.
if(CRPropa_FOUND AND NOT TARGET CRPropa::CRPropa)
	add_library(CRPropa::CRPropa UNKNOWN IMPORTED)
	set_target_properties(CRPropa::CRPropa PROPERTIES
		IMPORTED_LOCATION "${CRPropa_LIBRARY}"
		INTERFACE_INCLUDE_DIRECTORIES "${CRPropa_INCLUDE_DIR}"
		INTERFACE_LINK_DIRECTORIES "${CRPropa_LIBRARY_DIR}"
	)
endif()

if(CRPropa_FOUND AND NOT CRPropa_FIND_QUIETLY)
	message(STATUS "CRPropa include path: ${CRPropa_INCLUDE_DIR}")
	message(STATUS "CRPropa library: ${CRPropa_LIBRARY}")
	message(STATUS "CRPropa library path: ${CRPropa_LIBRARY_DIR}")
	message(STATUS "CRPropa SWIG interface file: ${CRPropa_SWIG_INTERFACE_FILE}")
	if(CRPropa_VERSION)
		message(STATUS "CRPropa version: ${CRPropa_VERSION}")
	endif()
endif()

mark_as_advanced(CRPropa_INCLUDE_DIR CRPropa_LIBRARY CRPropa_SWIG_PATH)
