# Define source files and add the executable.
project(cherenkov_library)
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
        monte_carlo.h
        data_analysis.cpp
        data_analysis.h)
add_library(cherenkov_library STATIC ${SOURCE_FILES})

# Link to external libraries.
include(ExternalLibraries.cmake)
link_boost(cherenkov_library)
link_root(cherenkov_library)
