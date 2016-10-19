# Set general flags and define the project
cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
project(cherenkov_library)

# Define source files and add the main executable
set(SOURCE_FILES
    accuracy_estimator.cpp
    common.cpp
    common.h
    data_containers.cpp
    data_containers.h
    geometric_objects.cpp
    geometric_objects.h
    reconstruction.cpp
    shower_simulator.cpp
    simulator.h)
add_library(cherenkov_library STATIC ${SOURCE_FILES})

# Link to the Boost package
list(APPEND CMAKE_PREFIX_PATH /usr/local/boost)
find_package(Boost COMPONENTS program_options REQUIRED)
target_link_libraries(cherenkov_library ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIR})

# Link to the ROOT package
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
target_link_libraries(cherenkov_library ${ROOT_LIBRARIES})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link to the FFTW library and header
target_link_libraries(cherenkov_library /usr/local/lib/libfftw3.a)
include_directories(/usr/local/include)

