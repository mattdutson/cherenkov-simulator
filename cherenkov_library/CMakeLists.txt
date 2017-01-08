# Set general flags and define the project
project(cherenkov_library)

# Define source files and add the main executable
set(SOURCE_FILES
        monte_carlo.cpp
        utility.cpp
        utility.h
        data_containers.cpp
        data_containers.h
        geometric_objects.cpp
        geometric_objects.h
        reconstructor.cpp
        simulator.cpp
        simulator.h
        reconstructor.h
        monte_carlo.h)
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
