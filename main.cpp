// main.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
// The entry point for the application

#include <fstream>
#include <iostream>
#include "cherenkov_library/simulator.h"

using namespace std;
using namespace cherenkov_simulator;

int main(int argc, const char* argv[])
{
    // Get the filename from the first command-line argument
    string filename = "config.xml";
    if (argc > 0)
    {
        filename = string(argv[0]);
    }

    // Try opening the specified file
    ifstream config_file = ifstream();
    try
    {
        config_file.open(filename);
    }
    catch (...)
    {
        cout << "The file " << filename << " could not be opened. Check the path." << endl;
        return -1;
    }

    // Try to construct the simulator using parameters in the configuration file
    Simulator sim = Simulator();
    try
    {
        sim.ParseFile(config_file);
    }
    catch (runtime_error e)
    {
        cout << e.what() << endl;
        return -1;
    }
    catch (...)
    {
        cout << "An unknown exception occured." << endl;
        return -1;
    }

    // Call methods on the simulator

    return 0;


}
