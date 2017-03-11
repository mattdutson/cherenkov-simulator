# Define source files and add the executable.
project(cherenkov_library)
set(SOURCE_FILES
        MonteCarlo.cpp
        Utility.cpp
        Utility.h
        DataStructures.cpp
        DataStructures.h
        Geometric.cpp
        Geometric.h
        Reconstructor.cpp
        Simulator.cpp
        Simulator.h
        Reconstructor.h
        MonteCarlo.h
        Analysis.cpp
        Analysis.h)
add_library(cherenkov_library STATIC ${SOURCE_FILES})

# Link to external libraries.
include(../external_lib/ExternalLib.cmake)
link_boost(cherenkov_library)
link_root(cherenkov_library)
