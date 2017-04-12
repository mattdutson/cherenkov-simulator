// Main.cpp
//
// Author: Matthew Dutson
//
// The entry point for the application.

#include "cherenkov_lib/MonteCarlo.h"

using boost::property_tree::ptree;

int main(int argc, const char* argv[])
{
    return cherenkov_simulator::MonteCarlo::Run(argc, argv);
}
