# GADGET_INCLUDE_DIR = path to gadget directory
# GADGET_LIBRARY = libGadget.a
# GADGET_FOUND = true if gadget is found

find_path(GADGET_INCLUDE_DIR gadget/SmoothParticle.h)
find_library(GADGET_LIBRARY Gadget)

set(GADGET_FOUND FALSE)
if(GADGET_INCLUDE_DIR AND GADGET_LIBRARY)
    set(GADGET_FOUND TRUE)
    MESSAGE(STATUS "Found Gadget")
else()
    MESSAGE(STATUS "Gadget NOT found")
endif()

MESSAGE(STATUS "  Include: ${GADGET_INCLUDE_DIR}")
MESSAGE(STATUS "  Library: ${GADGET_LIBRARY}")

mark_as_advanced(GADGET_INCLUDE_DIR GADGET_LIBRARY GADGET_FOUND)
