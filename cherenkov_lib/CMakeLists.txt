# Define source files and add the executable.
project(cherenkov_library)
set(SOURCE_FILES
        MonteCarlo.cpp
        Utility1.cpp
        Utility1.h
        DataStructures.cpp
        DataStructures.h
        Geometric.cpp
        Geometric.h
        Reconstructor1.cpp
        Simulator1.cpp
        Simulator1.h
        Reconstructor1.h
        MonteCarlo.h
        Analysis.cpp
        Analysis.h)
add_library(cherenkov_library STATIC ${SOURCE_FILES})

# Link to external libraries.
include(../external_lib/ExternalLib.cmake)
link_boost(cherenkov_library)
link_root(cherenkov_library)
