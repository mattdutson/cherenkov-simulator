// SimulatorTest.cpp
//
// Author: Matthew Dutson
//
// Tests of Simulator.h

#include <gtest/gtest.h>
#include <TFile.h>

#include "Helper.h"
#include "Utility.h"
#include "Simulator.h"
#include "MonteCarlo.h"
#include "Analysis.h"

using boost::property_tree::ptree;

namespace cherenkov_simulator
{
    class SimulatorTest : public testing::Test
    {

    protected:

        Simulator* simulator;
        MonteCarlo* monte_carlo;

        virtual void SetUp()
        {
            ptree config = Utility::ParseXMLFile("../../../Config.xml").get_child("config");
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
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        // Make some graphics and write them to a file.
        TFile file("StraightShowerSim.root", "RECREATE");
        Analysis::MakeTimeProfile(data).Write("straight_shower_graph");
        Analysis::MakePixlProfile(data).Write("straight_shower_map");
    }

    /*
     * Run a shower which goes sideways at some angle.
     */
    TEST_F(SimulatorTest, AngleShower)
    {
        // Construct and simulate a shower at 10km which is skewed at some angle.
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 0, -2), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        // Draw a map of impacts.
        TFile file("AngleShowerSim.root", "RECREATE");
        Analysis::MakeTimeProfile(data).Write("angle_shower_graph");
        Analysis::MakePixlProfile(data).Write("angle_shower_map");
    }

    /*
     * Run a shower which is representative of a typical shower.
     */
    TEST_F(SimulatorTest, TypicalShower)
    {
        // Construct and simulate a shower at 10km which moves slightly toward or away from the detector
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 1, -3), 1e6, -0.1, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        // Draw a map of impacts.
        TFile file("TypicalShowerSim.root", "RECREATE");
        Analysis::MakeTimeProfile(data).Write("typical_shower_graph");
        Analysis::MakePixlProfile(data).Write("typical_shower_map");
    }

    TEST_F(SimulatorTest, WholeProfile)
    {
        // Construct and simulate a shower at 30km.
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 3e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        // Make some graphics and write them to a file.
        TFile file("WholeProfileSim.root", "RECREATE");
        Analysis::MakeTimeProfile(data).Write("whole_shower_graph");
        Analysis::MakePixlProfile(data).Write("whole_shower_map");
    }
}