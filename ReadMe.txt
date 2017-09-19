Building the project requires several external libraries. Instructions for downloading and using them are given below.
Many of these instructions are specific to MacOS.

Google Test: https://github.com/google/googletest
The Google Test directory must be located at cherenkov_test/gtest. It is built as part of the test project.

ROOT:
There are two methods of installing ROOT: either directly as a tarball, or with the Homebrew package manager. In the
first method, ROOT can be downloaded from https://root.cern.ch/downloading-root with binaries and unpacked anywhere.
XCode and the command-line developer tools also need to be installed. XCode can be found on the Mac App Store, and the
command line tools can be installed with "xcode-select --install".

The disadvantage of simply downloading the ROOT tarball and extracting it is that it requires a full, bloated XCode
installation. Alternatively, ROOT can be installed with Homebrew (http://brew.sh). This only requires the XCode
command-line tools. To install, download Homebrew from the website and run "brew install root6". Note that this method
has resulted in some build issues in the past.

The thisroot.sh script sets environment variables and adds the ROOT executable to the path. The following line should be
added to the startup script (~/.bash_profile):
. ROOT_PATH/bin/thisroot.sh
ROOT_PATH should be replaced with the path to the ROOT installation (/usr/local/root).

Boost: http://www.boost.org/users/download/
Boost should be downloaded with precompiled binaries and can be placed anywhere on the system. The following line should
be added to the shell's startup script:
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:BOOST_PATH/stage/lib
This adds the Boost dynamic libraries to the search path. BOOST_PATH should be replaced by the path to the Boost
installation (/usr/local/boost).
