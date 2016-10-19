cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
project(cherenkov_tests)

set(SOURCE_FILES basic_tests.cpp)
add_executable(cherenkov_tests ${SOURCE_FILES})

include_directories(cherenkov_library)

add_subdirectory(gtest)

target_link_libraries(cherenkov_tests cherenkov_library)
target_link_libraries(cherenkov_tests gtest gtest_main)