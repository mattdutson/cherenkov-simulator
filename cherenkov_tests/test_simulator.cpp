//
// Created by Matthew Dutson on 1/19/17.
//

#include <gtest/gtest.h>
#include <vector>
#include <TFile.h>

#include "helper.h"
#include "utility.h"
#include "simulator.h"
#include "monte_carlo.h"
#include "data_analysis.h"

using boost::property_tree::ptree;
using std::vector;

using namespace cherenkov_library;

namespace cherenkov_tests
{
    class SimulatorTest : public testing::Test
    {

    protected:

        Simulator* simulator;
        MonteCarlo* monte_carlo;

        virtual void SetUp()
        {
            ptree config = ParseXMLFile("../../config.xml").get_child("config");
            simulator = new Simulator(config);
            monte_carlo = new MonteCarlo(config);
        }

        virtual void TearDown()
        {
            delete simulator;
            delete monte_carlo;
        }
    };

    /*
     * Attempt to simulate the motion of a typical shower to check whether profiles match expectations.
     */
    TEST_F(SimulatorTest, StraightShower)
    {
        // Construct and simulate a shower at 10 km. The axis coordinates are in the world frame.
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1000000, 0, 10e19);
        PhotonCount data = simulator->SimulateShower(shower);

        // Make some graphics and write them to a file.
        TFile file("../../cherenkov_tests/straight_shower_sim.root", "RECREATE");
        DataAnalysis::MakeProfileGraph(data).Write("straight_shower_graph");
        DataAnalysis::MakeSumMap(data).Write("straight_shower_map");
    }

    /*
     * Run a shower which goes sideways at some angle.
     */
    TEST_F(SimulatorTest, AngleShower)
    {
        // Construct and simulate a shower at 10km which is skewed at some angle.
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 0, -2), 1000000, 0, 10e19);
        PhotonCount data = simulator->SimulateShower(shower);

        // Draw a map of impacts.
        TFile file("../../cherenkov_tests/angle_shower_sim.root", "RECREATE");
        DataAnalysis::MakeProfileGraph(data).Write("angle_shower_graph");
        DataAnalysis::MakeSumMap(data).Write("angle_shower_map");
    }
}