# Find FFTW3 with single precision. Sets
# FFTW3F_FOUND = true if fftw3f is found
# FFTW3F_INCLUDE_DIR = fftw3.h
# FFTW3F_LIBRARY = libfftw3f.a .so

find_path(FFTW3F_INCLUDE_DIR fftw3.h)
find_library(FFTW3F_LIBRARY fftw3f)

set(FFTW3F_FOUND FALSE)
if(FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARY)
    set(FFTW3F_FOUND TRUE)
    MESSAGE("-- Found FFTW3 with single precision (FFTW3F):\n"
            "--   Include: ${FFTW3F_INCLUDE_DIR}\n"
            "--   Library: ${FFTW3F_LIBRARY}")
endif()

mark_as_advanced(FFTW3F_INCLUDE_DIR FFTW3F_LIBRARY FFTW3F_FOUND)