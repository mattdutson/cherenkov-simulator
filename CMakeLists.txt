cmake_minimum_required(VERSION 3.6)
project(cherenkov_simulator)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    accuracy_estimator.cpp
    common.cpp
    common.h
    data_containers.cpp
    data_containers.h
    geometric_objects.cpp
    geometric_objects.h
    main.cpp
    reconstruction.cpp
    shower_simulator.cpp
    simulator.h)

include_directories(/usr/local/boost /usr/local/root/include)

add_executable(cherenkov_simulator ${SOURCE_FILES})