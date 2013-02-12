cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(kiss STATIC
    src/logger
    src/path
    src/string
    src/uuid
)

SET_TARGET_PROPERTIES( kiss PROPERTIES COMPILE_FLAGS -fPIC)

# testing
if(ENABLE_TESTING)
    add_executable(test_uuid test/test_uuid.cpp)
    target_link_libraries(test_uuid kiss gtest gtest_main pthread)
    add_test(test_uuid test_uuid)
endif(ENABLE_TESTING)
