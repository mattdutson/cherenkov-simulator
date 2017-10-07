# Note that you'll want to set CMAKE_BUILD_TYPE to Release if you're looking for efficiency - in CLion this is done via
# a build configuration setting, and not directly in CMakeLists
cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Define source files and add the executable.
project(cherenkov_simulator)
set(SOURCE_FILES Main.cpp)
add_executable(cherenkov_simulator ${SOURCE_FILES})

# Link to cherenkov_lib and its dependencies.
include(ExternalLib.cmake)
link_boost(cherenkov_simulator)
link_root(cherenkov_simulator)
add_subdirectory(cherenkov_lib)
include_directories(cherenkov_lib)
target_link_libraries(cherenkov_simulator cherenkov_lib)

# Link to cherenkov_test and its dependencies.
add_subdirectory(cherenkov_test)

