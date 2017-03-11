# Find and link to the boost libraries.
function(link_boost library_name)
    list(APPEND CMAKE_PREFIX_PATH /usr/local/boost)
    find_package(Boost COMPONENTS program_options REQUIRED)
    target_link_libraries(${library_name} ${Boost_LIBRARIES})
    include_directories(${Boost_INCLUDE_DIR})
endfunction()

# Find and link to the ROOT libraries.
function(link_root library_name)
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
    find_package(ROOT REQUIRED COMPONENTS MathMore)
    include(${ROOT_USE_FILE})
    target_link_libraries(${library_name} ${ROOT_LIBRARIES})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
endfunction()

# Find and link to the FFTW libraries.
function(link_fftw library_name)
    target_link_libraries(${library_name} /usr/local/lib/libfftw3.a)
    include_directories(/usr/local/include)
endfunction()
