Building the project requires several external libraries. Instructions for downloading and using them are given below.

Google Test: https://github.com/google/googletest
The Google Test directory must be located at cherenkov_test/gtest. It is built as part of the test project.

ROOT: https://root.cern.ch/downloading-root
ROOT should be downloaded with precompiled binaries and can be placed anywhere on the system. The following line should
be added to the shell's startup script:
source ROOT_PATH/bin/thisroot.sh
The thisroot.sh script sets environment variables and adds the ROOT executable to the path. ROOT_PATH should be replaced
by the path to the ROOT installation (/usr/local/root).

Boost: http://www.boost.org/users/download/
Boost should be downloaded with precompiled binaries and can be placed anywhere on the system. The following line should
be added to the shell's startup script:
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:BOOST_PATH/stage/lib
This adds the Boost dynamic libraries to the search path. BOOST_PATH should be replaced by the path to the Boost
installation (/usr/local/boost).
