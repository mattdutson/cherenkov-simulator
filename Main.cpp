// Main.cpp
//
// Author: Matthew Dutson
//
// The entry point for the application.

#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "cherenkov_lib/Utility.h"

using namespace std;
using namespace cherenkov_lib;

using boost::property_tree::ptree;

int main(int argc, const char* argv[])
{
    // Get the filename from the first command-line argument
    string filename = "Config.xml";
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
        ptree config = cherenkov_lib::Utility::ParseXMLFile(filename).get_child("configuration");
    }
    catch (exception e)
    {
        cout << e.what() << endl;
        return -1;
    }
    return 0;
}
