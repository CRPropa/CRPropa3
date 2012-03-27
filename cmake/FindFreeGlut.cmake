#  try to find glut library and include files
#  Freeglut_INCLUDE_DIR, where to find GL/glut.h, etc.
#  Freeglut_LIBRARIES, the libraries to link against
#  Freeglut_FOUND, If false, do not try to use Freeglut.
# Also defined, but not for general use are:
#  Freeglut_glut_LIBRARY = the full path to the glut library.

#=============================================================================
# Copyright 2001-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

IF (WIN32)
  FIND_PATH( Freeglut_INCLUDE_DIR NAMES GL/freeglut.h 
    PATHS  ${Freeglut_ROOT_PATH}/include )
  FIND_LIBRARY( Freeglut_glut_LIBRARY NAMES glut glut32 freeglut freeglut_static
    PATHS
    ${OPENGL_LIBRARY_DIR}
    ${Freeglut_ROOT_PATH}/Release
    )
ELSE (WIN32)
  IF (APPLE)
    # These values for Apple could probably do with improvement.
    FIND_PATH( Freeglut_INCLUDE_DIR freeglut.h
      /System/Library/Frameworks/Freeglut.framework/Versions/A/Headers
      ${OPENGL_LIBRARY_DIR}
      )
    SET(Freeglut_glut_LIBRARY "-framework Freeglut" CACHE STRING "Freeglut library for OSX") 
    SET(Freeglut_cocoa_LIBRARY "-framework Cocoa" CACHE STRING "Cocoa framework for OSX")
  ELSE (APPLE)
    FIND_PATH( Freeglut_INCLUDE_DIR GL/freeglut.h
      /usr/include/GL
      /usr/openwin/share/include
      /usr/openwin/include
      /opt/graphics/OpenGL/include
      /opt/graphics/OpenGL/contrib/libglut
      )
  
    FIND_LIBRARY( Freeglut_glut_LIBRARY glut
      /usr/openwin/lib
      )
  ENDIF (APPLE)
ENDIF (WIN32)

SET( Freeglut_FOUND "NO" )
IF(Freeglut_INCLUDE_DIR)
  IF(Freeglut_glut_LIBRARY)
    SET( Freeglut_LIBRARIES
      ${Freeglut_glut_LIBRARY}
      ${Freeglut_cocoa_LIBRARY}
      )
    SET( Freeglut_FOUND "YES" )
    
  ENDIF(Freeglut_glut_LIBRARY)
ENDIF(Freeglut_INCLUDE_DIR)

MARK_AS_ADVANCED(
  Freeglut_INCLUDE_DIR
  Freeglut_glut_LIBRARY
  )