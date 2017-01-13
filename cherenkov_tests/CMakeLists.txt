# Define source files and add the executable.
project(cherenkov_tests)
set(SOURCE_FILES basic_tests.cpp)
add_executable(cherenkov_tests ${SOURCE_FILES})

# Link to cherenkov_library and its dependencies.
include_directories(../cherenkov_library)
include(../cherenkov_library/ExternalLibraries.cmake)
link_boost(cherenkov_tests)
link_root(cherenkov_tests)
target_link_libraries(cherenkov_tests cherenkov_library)

# Link to gtest.
add_subdirectory(gtest)
target_link_libraries(cherenkov_tests gtest gtest_main)

