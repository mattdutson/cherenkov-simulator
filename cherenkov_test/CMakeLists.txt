# Define source files and add the executable.
project(cherenkov_tests)
set(SOURCE_FILES MiscellaneousTest.cpp GeometricTest.cpp Helper1.h Helper1.cpp SimulatorTest.cpp DataStructuresTest.cpp ReconstructorTest.cpp)
add_executable(cherenkov_tests ${SOURCE_FILES})

# Link to cherenkov_lib and its dependencies.
include_directories(../cherenkov_lib)
include(../external_lib/ExternalLib.cmake)
link_boost(cherenkov_tests)
link_root(cherenkov_tests)
target_link_libraries(cherenkov_tests cherenkov_library)

# Link to gtest.
add_subdirectory(../external_lib/gtest)
target_link_libraries(cherenkov_tests ../external_lib/gtest gtest_main)
