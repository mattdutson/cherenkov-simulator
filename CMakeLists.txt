# Define general flags.
cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Define source files and add the executable.
project(cherenkov_simulator)
set(SOURCE_FILES main.cpp)
add_executable(cherenkov_simulator ${SOURCE_FILES})

# Add subdirectories.
include_directories(cherenkov_lib)
add_subdirectory(cherenkov_lib)
add_subdirectory(cherenkov_test)
include(external_lib/ExternalLib.cmake)
link_boost(cherenkov_simulator)
link_root(cherenkov_simulator)

# Link to the contents of cherenkov_lib.
target_link_libraries(cherenkov_simulator cherenkov_lib)
