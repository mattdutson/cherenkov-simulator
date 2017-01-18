// main.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
// The entry point for the application

#include <fstream>
#include <iostream>
#include "cherenkov_library/simulator.h"
#include "cherenkov_library/utility.h"

using namespace std;
using namespace cherenkov_library;

int main(int argc, const char* argv[])
{
    // Get the filename from the first command-line argument
    string filename = "config.xml";
    if (argc > 0)
    {
        filename = string(argv[0]);
    }
    try
    {
        ptree config = cherenkov_library::ParseXMLFile(filename).child("configuration");
    }
    catch (exception e)
    {
        cout << e.what() << endl;
        return -1;
    }

    // Call methods on the simulator

    return 0;
}
