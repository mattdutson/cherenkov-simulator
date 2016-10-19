cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
project(cherenkov_simulator)

set(SOURCE_FILES main.cpp)
add_executable(cherenkov_simulator ${SOURCE_FILES})

include_directories(cherenkov_library)

add_subdirectory(cherenkov_library)
add_subdirectory(cherenkov_tests)

target_link_libraries(cherenkov_simulator cherenkov_library)