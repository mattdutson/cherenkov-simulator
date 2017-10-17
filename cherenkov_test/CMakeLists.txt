set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Define source files and add the executable.
project(cherenkov_test)
set(SOURCE_FILES
        DataStructuresTest.cpp
        GeometricTest.cpp
        Helper.h
        Helper.cpp
        UtilityTest.cpp
        SampleEvents.cpp
        )
add_executable(cherenkov_test ${SOURCE_FILES})

# Link to cherenkov_lib and its dependencies.
include(../ExternalLib.cmake)
link_boost(cherenkov_test)
link_root(cherenkov_test)
include_directories(../cherenkov_lib)
target_link_libraries(cherenkov_test cherenkov_lib)

# Link to gtest.
add_subdirectory(gtest)
target_link_libraries(cherenkov_test gtest gtest_main)
