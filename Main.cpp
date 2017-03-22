// Main.cpp
//
// Author: Matthew Dutson
//
// The entry point for the application.

#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "cherenkov_lib/Utility.h"
#include "cherenkov_lib/MonteCarlo.h"

using namespace std;
using namespace cherenkov_simulator;

using boost::property_tree::ptree;

int main(int argc, const char* argv[])
{
    // Get the filename from the first command-line argument
    if (argc < 2)
    {
        cout << "Must specify output file as command line argument" << endl;
        return -1;
    }
    string out_file = argv[1];
    string config_file = "Config.xml";
    if (argc > 2) config_file = string(argv[2]);
    try
    {
        ptree config = cherenkov_simulator::Utility::ParseXMLFile(config_file).get_child("config");
        MonteCarlo monte_carlo = MonteCarlo(config);
        monte_carlo.PerformMonteCarlo(out_file);
        return 0;
    }
    catch (runtime_error err)
    {
        cout << err.what() << endl;
        return -1;
    }
}
