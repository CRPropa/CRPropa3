# - Finds ROOT instalation
# This module sets up ROOT information 
# It defines:
# ROOT_FOUND          If ROOT is found

find_program(ROOT_CONFIG_EXECUTABLE root-config
  PATHS $ENV{ROOTSYS}/bin)



if(NOT ROOT_CONFIG_EXECUTABLE)
  set(ROOT_FOUND FALSE)
  MESSAGE(STATUS "ROOT: NOT Found!")    
else()    
  set(ROOT_FOUND TRUE)
  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} "--cflags" OUTPUT_VARIABLE ROOT_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} "--libs" OUTPUT_VARIABLE ROOT_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "ROOT: Found!")
  message(STATUS "  CFlags:      " ${ROOT_CFLAGS})
  message(STATUS "  Libs:        " ${ROOT_LIBS})
endif()