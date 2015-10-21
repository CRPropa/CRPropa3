# SETUP PYTHON

# get default python inerpreter
FIND_PROGRAM( PYTHON_EXECUTABLE python REQUIRED
	PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath] [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.6\\InstallPath] [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.5\\InstallPath]
) 

SET(PYTHONINTERP_FOUND TRUE)

# find python include path
execute_process(
	COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(sysconfig.get_python_inc())"
	OUTPUT_VARIABLE PYTHON_INCLUDE_PATH
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

FIND_FILE(PYTHON_H_FOUND Python.h ${PYTHON_INCLUDE_PATH})
IF(NOT PYTHON_H_FOUND)
	MESSAGE(SEND_ERROR "Python.h not found")
ENDIF()

IF(APPLE)
	
	# apple changed linking to Python in OS X 10.9 (Mavericks)
	# until 10.8: use Python.framework as part of the SDK (-framework Python)
	# since 10.9: link to Python like any other UNIX

	# extract (minor) version number from SDK name (major version always 10 for OS X)
	STRING(REGEX REPLACE ".*MacOSX([0-9]+)[.]([0-9]+)[.]sdk" "\\2" OSX_SDK_MINOR_VERSION "${CMAKE_OSX_SYSROOT}" )
	# MESSAGE("Found OS X SDK minor version: ${OSX_SDK_MINOR_VERSION}")

	IF (PYTHON_LIBRARIES)
			SET(OSX_USE_PYTHON_FRAMEWORK "False")
			MESSAGE(STATUS "Using user provided Python library: " ${PYTHON_LIBRARIES} )

	ELSE(PYTHON_LIBRARIES)
		IF(OSX_SDK_MINOR_VERSION GREATER 8)
			SET(OSX_USE_PYTHON_FRAMEWORK "False")
			MESSAGE(STATUS "Running on Mac OS X >= 10.9: Linking to Python in UNIX default way")
		ELSE(OSX_SDK_MINOR_VERSION GREATER 8)
			MESSAGE(STATUS "Running on Mac OS X < 10.9: Linking to Python as framework")
			SET(OSX_USE_PYTHON_FRAMEWORK "True")

			INCLUDE(CMakeFindFrameworks)
			# Search for the python framework on Apple.
			MESSAGE(INFO "Looking for python framework as on apple system" )
			CMAKE_FIND_FRAMEWORKS(Python)
			SET (PYTHON_LIBRARIES "-framework Python" CACHE FILEPATH "Python Framework" FORCE)
		ENDIF(OSX_SDK_MINOR_VERSION GREATER 8)
	ENDIF(PYTHON_LIBRARIES)

ENDIF(APPLE)

IF(MSVC)
	execute_process(
		COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; import os; prefix= sysconfig.get_config_var('prefix'); ver = sysconfig.get_python_version().replace('.', ''); lib = os.path.join(prefix,'libs\\python'+ver+'.lib'); sys.stdout.write(lib)"
		OUTPUT_VARIABLE PYTHON_LIBRARIES
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
ENDIF(MSVC)

IF (MINGW)
	execute_process(
		COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; import os; prefix= sysconfig.get_config_var('prefix'); ver = sysconfig.get_python_version().replace('.', ''); lib = os.path.join(prefix,'libs\\libpython'+ver+'.a'); sys.stdout.write(lib)"
		OUTPUT_VARIABLE PYTHON_LIBRARIES
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
ENDIF(MINGW)

IF(NOT OSX_USE_PYTHON_FRAMEWORK AND NOT PYTHON_LIBRARIES AND NOT MSVC AND NOT MINGW)
	execute_process(
		COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; import os; libname = sysconfig.get_config_var('LDLIBRARY'); libdir = sysconfig.get_config_var('LIBPL'); lib = os.path.join(libdir,libname); out = lib if os.path.exists(lib) else os.path.join(libdir, sysconfig.get_config_var('LIBRARY')); sys.stdout.write(out);"
		OUTPUT_VARIABLE PYTHON_LIBRARIES
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
ENDIF(NOT OSX_USE_PYTHON_FRAMEWORK AND NOT PYTHON_LIBRARIES AND NOT MSVC AND NOT MINGW)

#find the site package destinaton 
execute_process(
	COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
	OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
	COMMAND ${PYTHON_EXECUTABLE} -c "import sys; sys.stdout.write(str(sys.version_info[0]))"
	OUTPUT_VARIABLE PYTHON_MAJOR_VERSION
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
	COMMAND ${PYTHON_EXECUTABLE} -c "import sys; sys.stdout.write(str(sys.version_info[1]))"
	OUTPUT_VARIABLE PYTHON_MINOR_VERSION
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

SET(PYTHON_VERSION "${PYTHON_MAJOR_VERSION}${PYTHON_MINOR_VERSION}")
SET(PYTHON_DOT_VERSION "${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION}")

MESSAGE(STATUS "Python: Found!")
MESSAGE(STATUS "  Version:     " ${PYTHON_DOT_VERSION} "/" ${PYTHON_VERSION})
MESSAGE(STATUS "  Executeable: " ${PYTHON_EXECUTABLE})
MESSAGE(STATUS "  Include:     " ${PYTHON_INCLUDE_PATH})
MESSAGE(STATUS "  Library:     " ${PYTHON_LIBRARIES})
MESSAGE(STATUS "  Site-package directory: " ${PYTHON_SITE_PACKAGES})
