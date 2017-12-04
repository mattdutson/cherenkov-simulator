# Cherenkov Simulator

## Purpose and Background
This code was written by Matthew Dutson at the University of Utah under the direction of Prof. Douglas Bergman. It is part of the University of Utah's ongoing research into ultra high energy cosmic rays. More information about cosmic ray research and the Telescope Array Project can be found [here](http://www.telescopearray.org).

This simulation investigates a possible improvement to monocular air fluorescence reconstruction. In the proposed method, the reflection point of Cherenkov radiation is used to reduce the number of free parameters in a time profile fit. The details of the method and the theory behind this code are given [here](http://bit.ly/2jcnh24).

The core of this application is a Monte Carlo simulation. The Monte Carlo randomly generates cosmic ray cascades within a certain parameter space. It simulates the propagation of this cascade and its detection by a fluorescence telescope. The improved reconstruction method is used and compared against the traditional method. Over thousands of iteration, this gives a good comparison of performance.

## Required Libraries
Builds have only been attempted on macOS 10.12 and 10.13. Adjustments may be required on other platforms.

### Google Test
Google Test is used both for unit testing and for simulating individual test showers. [Download](https://github.com/google/googletest) the Google Test source and place it at `cherenkov_test/gtest`. It is built simultaneously with the `cherenkov_test` project.

### Boost
Boost doesn't require any binaries. Headers should be [downloaded](http://www.boost.org/users/download/) and extracted to `external_lib/boost`.

### ROOT
ROOT should be [downloaded](https://root.cern.ch/downloading-root) with OS-specific binaries. On macOS, it requires that XCode, the XCode command-line tools, and XQuartz first be installed. The unzipped ROOT directory should be placed at `external_lib/root`. Locate the environment script `root/bin/thisroot.sh` and set it to run upon shell startup (probably by modifying `~/.bash_profile`). The script sets environment variables and adds the ROOT executable to the path. If this fails to correctly reference dynamic libraries, append `root/lib` to your `DYLD_LIBRARY_PATH`.
