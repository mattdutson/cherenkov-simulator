project(cherenkov_tests)

set(SOURCE_FILES basic_tests.cpp)
add_executable(cherenkov_tests ${SOURCE_FILES})

include_directories(cherenkov_library)

add_subdirectory(gtest)

target_link_libraries(cherenkov_tests gtest gtest_main)
target_link_libraries(cherenkov_tests cherenkov_library)
