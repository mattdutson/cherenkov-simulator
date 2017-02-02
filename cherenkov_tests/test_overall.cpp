//
// Created by Matthew Dutson on 1/19/17.
//

#include <gtest/gtest.h>
#include <vector>
#include <TFile.h>

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
        Simulator simulator = Simulator(config);
        MonteCarlo monte_carlo = MonteCarlo(config);

        // Construct and simulate a shower at 10 km. The axis coordinates are in the world frame.
        Shower shower = monte_carlo.GenerateShower(TVector3(0, 0, -1), 1000000, 0, 10e19);
        PhotonCount data = simulator.SimulateShower(shower);

        // Make some graphics and write them to a file.
        TGraph graph = MakeProfileGraph(data);
        TH2I histo = MakeSumMap(data);
        TFile file("../../cherenkov_tests/sample_shower.root", "RECREATE");
        graph.Write("straight_down");
        histo.Write("shower_map");
    }

    /*
     * Run a shower which goes sideways at some angle.
     */
    TEST(sample_shower, angle_motion)
    {
        // Parse the XML file and give it to the simulator and monte carlo.
        ptree config = ParseXMLFile("../../config.xml").get_child("config");
        Simulator simulator = Simulator(config);
        MonteCarlo monte_carlo = MonteCarlo(config);

        // Construct and simulate a shower at 10km which is skewed at some angle.
        Shower shower = monte_carlo.GenerateShower(TVector3(1, 0, -2), 1000000, 0, 10e19);
        PhotonCount data = simulator.SimulateShower(shower);

        // Draw a map of impacts.
        TGraph graph = MakeProfileGraph(data);
        TH2I histo = MakeSumMap(data);
        TFile file("../../cherenkov_tests/angle_shower.root", "RECREATE");
        graph.Write("angle_shower_graph");
        histo.Write("angle_shower_histo");
    }
}