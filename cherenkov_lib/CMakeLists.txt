# Define source files and add the executable.
project(cherenkov_lib)
set(SOURCE_FILES
        Analysis.cpp
        Analysis.h
        DataStructures.cpp
        DataStructures.h
        Geometric.cpp
        Geometric.h
        MonteCarlo.cpp
        MonteCarlo.h
        Reconstructor.cpp
        Reconstructor.h
        Simulator.cpp
        Simulator.h
        Utility.cpp
        Utility.h)
add_library(cherenkov_lib STATIC ${SOURCE_FILES})

# Link to external libraries.
include(../external_lib/ExternalLib.cmake)
link_boost(cherenkov_lib)
link_root(cherenkov_lib)
