//
// Created by Matthew Dutson on 1/19/17.
//

#include <gtest/gtest.h>
#include <vector>

#include "test_helper.h"
#include "utility.h"
#include "simulator.h"
#include "monte_carlo.h"
#include "data_analysis.h"

using boost::property_tree::ptree;
using std::vector;

using namespace cherenkov_library;

namespace cherenkov_tests
{
    /*
     * Attempt to simulate the motion of a typical shower to check whether profiles match expectations.
     */
    TEST(sample_shower, straight_down)
    {
        // Parse the XML file and give it to the simulator and monte carlo.
        ptree config = ParseXMLFile("../../config.xml").get_child("config");
        Simulator simulator = Simulator();
        simulator.ParseFile(config);
        MonteCarlo monte_carlo = MonteCarlo();
        monte_carlo.ParseFile(config);

        // Construct and simulate a shower. The axis coordinates are in the world frame.
        Shower shower = monte_carlo.GenerateShower(TVector3(0, 0, -1), 100000, 0);
        PhotonCount data = simulator.SimulateShower(shower);

        // Output data to a CSV file.
        vector<vector<double>> profile = SuperimposeTimes(data);
        WriteCSV(profile, {"Time", "Count"}, "../../shower_output.csv");
    }
}