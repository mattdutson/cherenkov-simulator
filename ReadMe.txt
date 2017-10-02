Building the project requires several external libraries. Instructions for downloading and using them are given below.

Google Test: https://github.com/google/googletest
The Google Test directory must be located at cherenkov_test/gtest. It is built as part of the test project.

ROOT: https://root.cern.ch/downloading-root
ROOT should be downloaded with OS-specific binaries. On Mac OS, it requires that XCode, the XCode command-line tools,
and XQuartz first be installed. The unzipped ROOT directory should be placed at PROJECT_ROOT/external_lib/root. The
ExternalLib.cmake file searches the external_lib directory when linking. The thisroot.sh script sets environment
variables and adds the ROOT executable to the path. The following line should be added to the shell's startup script:
. PROJECT_ROOT/external_lib/root/bin/thisroot.sh

Boost: http://www.boost.org/users/download/
The Boost source should be downloaded and unzipped into PROJECT_ROOT/external_lib/boost.
