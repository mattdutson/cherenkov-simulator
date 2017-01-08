// common.h
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
// Contains a configuration manager and globally useful methods

#ifndef Common_h
#define Common_h

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"

namespace cherenkov_simulator
{
    /*
     * Reads the string and converts it to a TVector3 object.
     */
    TVector3 ToVector(std::string s);

    /*
     * Generates a randomly rotated vector perpendicular to the input. If the input vector is zero, (1, 0, 0) is
     * returned.
     */
    TVector3 RandomPerpendicularVector(TVector3 vec, TRandom3 rng);
}


#endif
