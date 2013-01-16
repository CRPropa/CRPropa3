# - Finds ROOT instalation
# This module sets up ROOT information 
# It defines:
# ROOT_FOUND          If the ROOT is found

find_program(ROOT_CONFIG_EXECUTABLE root-config
  PATHS $ENV{ROOTSYS}/bin)

if(NOT ROOT_CONFIG_EXECUTABLE)
  set(ROOT_FOUND FALSE)
  MESSAGE(STATUS "ROOT: NOT Found!")    
else()    
  set(ROOT_FOUND TRUE)
  MESSAGE(STATUS "ROOT: Found!")
endif()