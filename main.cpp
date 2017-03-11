// main.cpp
// cherenkov_lib
//
// Created by Matthew Dutson on 9/8/16.
//
// The entry point for the application

#include <fstream>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "cherenkov_library/simulator.h"
#include "cherenkov_library/monte_carlo.h"
#include "cherenkov_library/reconstructor.h"
#include "cherenkov_library/utility.h"

using namespace std;
using namespace cherenkov_library;

using boost::property_tree::ptree;

int main(int argc, const char* argv[])
{
    // Get the filename from the first command-line argument
    string filename = "Config1.xml";
    if (argc > 0)
    {
        filename = string(argv[0]);
    }
    else
    {
        return 0;
    }
    try
    {
        ptree config = cherenkov_library::ParseXMLFile(filename).get_child("configuration");
    }
    catch (exception e)
    {
        cout << e.what() << endl;
        return -1;
    }
    return 0;
}
