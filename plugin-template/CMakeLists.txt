cmake_minimum_required(VERSION 3.14)

# choose your plugin name
set(pluginName "myPlugin")

# set the version to whatever you want
project(${pluginName} 
	VERSION 1.0.0 
	LANGUAGES Fortran C CXX
)



###############################################################################
#                   Compiler- and system-specific options                     #
###############################################################################

# C++ should be at least 17
set(CMAKE_CXX_STANDARD 17)

# set default build-type to release to enable performance improvements
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

# set other needed compiler flags:
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")

# enable extra warnings on debug builds
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wextra")

# find cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# compiler-specific flags for Linux and macOS
if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
	message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)

# RPATH handling (Linux and macOS):
#   The plugin (and its Python module) must find libcrpropa, whose install name is @rpath/libcrpropa.
#   They must also find their own libraries at runtime.
#   Here we:
#   * keep the linker's build rpath so binaries run from the build tree, and
#   * add both the plugin's own install lib dir and the directories of all linked libraries (CRPropa included) to the install rpath.
set(CMAKE_MACOSX_RPATH On)
set(CMAKE_SKIP_BUILD_RPATH False)
set(CMAKE_BUILD_WITH_INSTALL_RPATH False)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

# display compiler information
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER} (version ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "C++ compilation flags: ${CMAKE_CXX_FLAGS}")

# general options
option(ENABLE_PYTHON "Create Python bindings with SWIG" On)
option(USE_OWN_CRPROPA "Use builtin CRPropa instead of in CRPROPA_PATH installed CRPropa" Off)
option(ENABLE_TESTING "Build tests and enable test target" On)
option(ENABLE_COVERAGE "Add code for coverage analysis." Off)
option(ENABLE_PYBIND11_STUBGEN "Generate Pybind11 stubs" Off)


# This token means "relative to the loaded binary" on each platform.
# It should work for unix-based platforms (Linux, MacOS).
if(APPLE)
	set(originToken "@loader_path")
else()
	set(originToken "$ORIGIN")
endif(APPLE)

# search next to the binary and in the plugin's own install lib directory
set(CMAKE_INSTALL_RPATH "${originToken};${CMAKE_INSTALL_PREFIX}/lib")
message(STATUS "Install RPATH: ${CMAKE_INSTALL_RPATH} (+ linked library dirs)")


###############################################################################
#                   Compilation and installation of plugin                    #
###############################################################################

# finding Python
if(ENABLE_PYTHON)
	set(Python_FIND_VIRTUALENV FIRST)
	set(Python_FIND_STRATEGY LOCATION)
	find_package(Python 3 REQUIRED COMPONENTS Interpreter Development.Module NumPy)
	
	if(Python_FOUND)
		message(STATUS "Python version: ${Python_VERSION}")
		message(STATUS "Python executable: ${Python_EXECUTABLE}")
		message(STATUS "Python headers: ${Python_INCLUDE_DIRS}")
		message(STATUS "Python library: ${Python_LIBRARY}")
		message(STATUS "Python site-packages: ${Python_SITE_PACKAGES}")
		add_definitions(-D${pluginName}_HAVE_PYTHON)
	endif(Python_FOUND)

	# find SWIG
	find_package(SWIG REQUIRED)
	include(UseSWIG) 

	# only support versions of SWIG higher than 4.0
	if(SWIG_VERSION VERSION_LESS 4.0)
		message(STATUS "SWIG version too old (${SWIG_VERSION} < 4.0). Code will not be compile with Python support.")
		set(SWIG_FOUND False)
	endif(SWIG_VERSION VERSION_LESS 4.0)	

	if(SWIG_FOUND)
		message(STATUS "SWIG installation directory: ${SWIG_DIR}")
		message(STATUS "SWIG executable: ${SWIG_EXECUTABLE}")
		list(APPEND ${pluginName}_SWIG_DEFINES -D${pluginName}_HAVE_PYTHON)
	else()
		message(STATUS "SWIG: NOT FOUND (cannot build with Python bindings)")
	endif(SWIG_FOUND)
	
	# find NumPy
	if(Python_NumPy_FOUND)
		include_directories(${Python_NumPy_INCLUDE_DIRS})
		message(STATUS "NumPy headers: ${Python_NumPy_INCLUDE_DIRS} (version ${Python_NumPy_VERSION})")
	else()
		message(STATUS "NumPy: NOT FOUND.")
		message(STATUS "The code might work just fine without NumPy, but some features in Python might not be straightfoward to use.")
	endif(Python_NumPy_FOUND)

	#  use Python_INSTALL_PACKAGE_DIR if provided; otherwise, install in Python_SITELIB
	set(Python_INSTALL_PACKAGE_DIR "${Python_SITELIB}" CACHE PATH "folder in which the python package is installed")
	message(STATUS "  package install directory: ${Python_INSTALL_PACKAGE_DIR}")

endif(ENABLE_PYTHON)

# support coverage
if(ENABLE_COVERAGE)
	if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
		set(COVERAGE_LIBS "-coverage -lgcov")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
	endif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
endif(ENABLE_COVERAGE)


# add other extra libraries by appending to these variables (see below on the example of CRPropa):
set(${pluginName}_EXTRA_INCLUDES)
set(${pluginName}_EXTRA_LIBRARIES)
set(${pluginName}_EXTRA_SOURCES)


# build crpropa and add resulting libcrpropa.so as library
if(USE_OWN_CRPROPA)
	if(NOT crpropa_POPULATED)
		include(FetchContent)
		
		# pulls repository specified and can even be pinned on specific commit
		FetchContent_Declare(
			crpropa
			GIT_REPOSITORY "https://github.com/CRPropa/CRPropa3.git"
			# # The GIT_TAG option enables you to specify a commit hash or specific tag, see:
			# #   https://cmake.org/cmake/help/latest/module/ExternalProject.html#git
			# GIT_TAG f5cf064cf60b8a4771deab6d963b01d44df13c3f
			GIT_TAG master
		)

		# contains .cmake files necessary for finding modules
		list(APPEND CMAKE_MODULE_PATH "${CRPropa_SOURCE_DIR}/cmake")
		
		# # may be useful to set some options for CRPropa before building it, e.g.:
		# set(CMAKE_ARGS per set before add_subdirectory (FORCE option overwrites existing entries)
		# set(ENABLE_GIT OFF CACHE BOOL "")
		# set(OMP_SCHEDULE "dynamic" CACHE STRING "")
		# set(FAST_WAVES ON CACHE BOOL "")
		# set(SIMD_EXTENSIONS "avx+fma" CACHE STRING "")
		
		FetchContent_MakeAvailable(crpropa)
	endif(NOT crpropa_POPULATED)

	set(CRPropa_SOURCE_DIR "${crpropa_SOURCE_DIR}")
	set(CRPropa_BINARY_DIR "${crpropa_BINARY_DIR}")

	list(APPEND ${pluginName}_EXTRA_INCLUDES ${CRPropa_SOURCE_DIR}/include)
	list(APPEND ${pluginName}_EXTRA_INCLUDES ${CRPropa_BINARY_DIR}/include)
	list(APPEND ${pluginName}_EXTRA_INCLUDES ${CRPROPA_EXTRA_INCLUDES})
	list(APPEND ${pluginName}_EXTRA_LIBRARIES crpropa)
	set(CRPropa_SWIG_INTERFACE_FILES ${CRPropa_SOURCE_DIR}/python)
else()
	# Find an installed CRPropa using cmake/FindCRPropa.cmake (found via CMAKE_MODULE_PATH set above). 
	# Every location can be overridden manually:
	# - `CRPropa_INSTALL_PREFIX`: the prefix where CRPropa is installed (e.g. /user/CRPropa3/build)
	# - `CRPropa_INCLUDE_DIR`: the directory containing CRPropa headers
	# - `CRPropa_LIBRARY`: the CRPropa library
	# - `CRPropa_SWIG_PATH`: the path to CRPropa SWIG interface files
	find_package(CRPropa REQUIRED)

	# append header to includes, library to libraries and remember location of sw ig_interface files
	list(APPEND ${pluginName}_EXTRA_INCLUDES ${CRPropa_INCLUDE_DIR})
	list(APPEND ${pluginName}_EXTRA_LIBRARIES ${CRPropa_LIBRARY})
	set(CRPropa_SWIG_INTERFACE_FILES ${CRPropa_SWIG_PATH})
endif(USE_OWN_CRPROPA)


# Add the source files located in ${CMAKE_CURRENT_SOURCE_DIR}/src
add_library(${pluginName} SHARED
	${CMAKE_CURRENT_SOURCE_DIR}/src/plugin.cpp
	${${pluginName}_EXTRA_SOURCES}
)

# the plugin's own headers plus any extra includes (e.g. CRPropa's CMakeLists.txt)
target_include_directories(${pluginName} PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${${pluginName}_EXTRA_INCLUDES}
)

# link crpropa and other needed libraries
target_link_libraries(${pluginName} ${${pluginName}_EXTRA_LIBRARIES})

# Python interface, this part should work basically fully automated.
if(ENABLE_PYTHON)
	if(Python_FOUND)
		# use SWIG's builtin (maybe should be abandoned if everything will be tested and working with -builtin)
		option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" On)

		if(ENABLE_SWIG_BUILTIN)
			set(BUILTIN "-builtin")
			list(APPEND ${pluginName}_SWIG_DEFINES "-builtin")
		else()
			set(BUILTIN "")
		endif(ENABLE_SWIG_BUILTIN)

		file(GLOB_RECURSE ${pluginName}_SWIG_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/python/*.i)
		set(${pluginName}_PYTHON_PACKAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/${pluginName}")
		file(MAKE_DIRECTORY "${${pluginName}_PYTHON_PACKAGE_DIR}")
		file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/python/${pluginName}/" DESTINATION "${${pluginName}_PYTHON_PACKAGE_DIR}")

		set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/python/${pluginName}.i
			PROPERTIES
			CPLUSPLUS ON
			SWIG_MODULE_NAME "${pluginName}"
			SWIG_FLAGS "${${pluginName}_SWIG_DEFINES};-dirprot;${CMAKE_SWIG_FLAGS}"
			OUTPUT_DIR "${${pluginName}_PYTHON_PACKAGE_DIR}"
			INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/include;${CRPropa_SWIG_INTERFACE_FILES};${${pluginName}_EXTRA_INCLUDES};"
		)
		
		swig_add_library(${pluginName}-swig TYPE SHARED LANGUAGE python
			SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/python/${pluginName}.i
		)
		
		target_include_directories(${pluginName}-swig PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}/include
			${${pluginName}_EXTRA_INCLUDES}
			${Python_INCLUDE_DIRS}
			${Python_NumPy_INCLUDE_DIRS}
			${CRPropa_SWIG_INTERFACE_FILES}
		)
		
		target_link_libraries(${pluginName}-swig ${pluginName} ${Python_LIBRARIES} ${Python_LIBRARY} Python::Module)

		
		set_target_properties(${pluginName}-swig PROPERTIES
			PREFIX ""
			OUTPUT_NAME "_${pluginName}"
			LIBRARY_OUTPUT_DIRECTORY "${${pluginName}_PYTHON_PACKAGE_DIR}"
			COMPILE_FLAGS  "-w"
			SUFFIX ".so"
		)

		if(APPLE)
			set_target_properties(${pluginName}-swig PROPERTIES INSTALL_RPATH "@loader_path;${CMAKE_INSTALL_FULL_LIBDIR}")
			set_property(TARGET ${pluginName}-swig APPEND PROPERTY LINK_FLAGS "-flat_namespace -undefined suppress")
		else()
			set_target_properties(${pluginName}-swig PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
		endif(APPLE)

	endif(Python_FOUND)

endif(ENABLE_PYTHON)


###############################################################################
#                                 Testing                                     #
###############################################################################

# the testing part, where tests for your code should be added
if(ENABLE_TESTING)
	enable_testing()
	set(${pluginName}_TESTS_PATH ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "")  
	mark_as_advanced(${pluginName}_TESTS_PATH)


	# temporary workaround for newer clang versions due to its handling of unwinding
	#   see: https://github.com/google/googletest/issues/3062
	if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
		if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
			add_definitions(-D${pluginName}_TESTS_SKIP_EXCEPTIONS)
		endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
	endif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

	if(APPLE)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
	endif(APPLE)

	# Reuse CRPropa's own gtest if at all possible; only if that fails, look for the system's gtest (or download it)
	set(${pluginName}_GTEST_LIBRARIES)
	set(${pluginName}_GTEST_INCLUDE_DIRS)

	# install CRPropa together with plugin
	if(USE_OWN_CRPROPA)
		set(${pluginName}_GTEST_INCLUDE_DIRS "${CRPropa_SOURCE_DIR}/libs/gtest/include")
		set(${pluginName}_GTEST_LIBRARIES gtest gtest_main)

	else()
		# check if CRPropa's bundled gtest is available (it should be, if CRPropa was built with ENABLE_TESTING, for versions of CRPropa >= 3.3) 
		# Resolve GTest from CRPropa's bundle first, else fall back to system.
		find_library(gtestCRPropa     NAMES gtest      PATHS "${CRPropa_LIBRARY_DIR}" NO_DEFAULT_PATH)
		find_library(gtestCRPropaMain NAMES gtest_main PATHS "${CRPropa_LIBRARY_DIR}" NO_DEFAULT_PATH)
		if(IS_DIRECTORY "${CRPropa_INCLUDE_DIR}/gtest" AND gtestCRPropa AND gtestCRPropaMain)
			set(${pluginName}_GTEST_INCLUDE_DIRS "${CRPropa_INCLUDE_DIR}")
			set(${pluginName}_GTEST_LIBRARIES    "${gtestCRPropa}" "${gtestCRPropaMain}")
			message(STATUS "GTest: using CRPropa's bundled version, (${gtestCRPropa}, ${gtestCRPropaMain})")
		# otherwise try to find a system-installed gtest package (or download it if not found)
		else()
			find_package(GTest QUIET)
			if(TARGET GTest::gtest)
				set(${pluginName}_GTEST_LIBRARIES GTest::gtest)
				message(STATUS "GTest: using system's GTest package")
			elseif(TARGET GTest::GTest)
				set(${pluginName}_GTEST_LIBRARIES GTest::GTest)
				message(STATUS "GTest: using system's legacy GTest package.")
			else()
				include(FetchContent)
				FetchContent_Declare(
					googletest
					GIT_REPOSITORY https://github.com/google/googletest.git
					GIT_TAG v1.14.0 # v1.17.0 fails on some systems
				)
				FetchContent_MakeAvailable(googletest)
				set(${pluginName}_GTEST_LIBRARIES gtest)
				list(APPEND ${pluginName}_GTEST_INCLUDE_DIRS "${gtest_SOURCE_DIR}/include" "${gtest_SOURCE_DIR}")
				message(STATUS "Downloading GTests from Github.")
			endif()
		endif(IS_DIRECTORY "${CRPropa_INCLUDE_DIR}/gtest" AND gtestCRPropa AND gtestCRPropaMain)
	endif(USE_OWN_CRPROPA)

	# function to create C++ tests, which are compiled and linked against gtest and the plugin
	function(addTest testName)
		add_executable(${testName} ${CMAKE_CURRENT_SOURCE_DIR}/test/${testName}.cpp)
		target_include_directories(${testName} PRIVATE ${${pluginName}_GTEST_INCLUDE_DIRS})
		target_link_libraries(${testName} PRIVATE ${pluginName} ${${pluginName}_GTEST_LIBRARIES} pthread ${COVERAGE_LIBS})
		add_test(NAME ${testName}Cpp COMMAND ${${pluginName}_TESTS_PATH}/${testName})
	endfunction()

	# first add the cpp tests
	addTest(testPlugin testPlugin)



	if(ENABLE_PYTHON AND Python_FOUND)

		file(
			COPY ${CMAKE_CURRENT_SOURCE_DIR}/test/
			DESTINATION  ${CMAKE_CURRENT_BINARY_DIR}
			FILES_MATCHING PATTERN "*.py"
		)

		if(USE_OWN_CRPROPA)
			set(CRPropa_PYTHON_TEST_PACKAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/crpropa")
			file(MAKE_DIRECTORY "${CRPropa_PYTHON_TEST_PACKAGE_DIR}")
			file(WRITE "${CRPropa_PYTHON_TEST_PACKAGE_DIR}/__init__.py" "from .crpropa import *\n")
			add_test(
				NAME moveCRPropaPlugin
				COMMAND ${CMAKE_COMMAND} -E copy
					"${crpropa_BINARY_DIR}/crpropa.py"
					"$<TARGET_FILE:crpropa-swig>"
					"$<TARGET_FILE:crpropa>"
					"${CRPropa_PYTHON_TEST_PACKAGE_DIR}"
			)
		endif(USE_OWN_CRPROPA)

		function(addTestPy testName)
			configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/${testName}.py ${CMAKE_CURRENT_BINARY_DIR}/${testName}.py COPYONLY)
			add_test(NAME ${testName}Py COMMAND ${Python_EXECUTABLE} ${${pluginName}_TESTS_PATH}/${testName}.py)
			if(USE_OWN_CRPROPA)
				set_tests_properties(${testName}Py PROPERTIES DEPENDS moveCRPropaPlugin)
			endif(USE_OWN_CRPROPA)
		endfunction()

		# add python tests here, just copy paste the add_test statements and change test name and path
		addTestPy(testPlugin)


	endif(ENABLE_PYTHON AND Python_FOUND)

endif(ENABLE_TESTING)



###############################################################################
#                                Installation                                 #
###############################################################################

# specify installation instructions
install(TARGETS ${pluginName} DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION "${CMAKE_INSTALL_PREFIX}/include" FILES_MATCHING PATTERN "*.h")
if(ENABLE_PYTHON AND Python_FOUND)
	install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/${pluginName}" DESTINATION "${Python_INSTALL_PACKAGE_DIR}")
	install(FILES "${${pluginName}_PYTHON_PACKAGE_DIR}/${pluginName}.py" DESTINATION "${Python_INSTALL_PACKAGE_DIR}/${pluginName}")
	install(TARGETS ${pluginName}-swig LIBRARY DESTINATION "${Python_INSTALL_PACKAGE_DIR}/${pluginName}")
	install(FILES ${${pluginName}_SWIG_INPUTS} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/${pluginName}/swig_interface)
endif(ENABLE_PYTHON AND Python_FOUND)


if(ENABLE_PYBIND11_STUBGEN AND ENABLE_PYTHON AND Python_FOUND)
	find_program(PYBIND11_STUBGEN_EXECUTABLE NAMES "pybind11-stubgen")
	if(PYBIND11_STUBGEN_EXECUTABLE)
		message(STATUS "Found pybind11-stubgen at ${PYBIND11_STUBGEN_EXECUTABLE}")
		install(CODE "execute_process(COMMAND ${PYBIND11_STUBGEN_EXECUTABLE} --ignore-all-errors -o ${Python_SITELIB} ${pluginName})")
		install(CODE "execute_process(COMMAND echo -- Generated python stubs over pybind11-stubgen for ${pluginName})")
	endif(PYBIND11_STUBGEN_EXECUTABLE)
endif(ENABLE_PYBIND11_STUBGEN AND ENABLE_PYTHON AND Python_FOUND)
