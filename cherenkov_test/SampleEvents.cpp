// SampleEvents.cpp
//
// Author: Matthew Dutson
//
// Can run several different sample showers.

#include <gtest/gtest.h>
#include <TFile.h>

#include "Reconstructor.h"
#include "Simulator.h"
#include "MonteCarlo.h"
#include "Analysis.h"

using namespace boost::property_tree;
using namespace std;
using namespace TMath;

namespace cherenkov_simulator
{
    class SampleEvents : public testing::Test
    {
    protected:

        MonteCarlo* monte_carlo;

        virtual void SetUp()
        {
            ptree config = Utility::ParseXMLFile("../Config.xml").get_child("config");
            monte_carlo = new MonteCarlo(config);
        }

        virtual void TearDown()
        {
            delete monte_carlo;
        }

        Plane FriendGroundPlane()
        {
            return monte_carlo->simulator.GroundPlane();
        }
    };

    /*
     * Simulate and attempt reconstruction on a straight vertical shower at 10km.
     */
    TEST_F(SampleEvents, StraightShower)
    {
        TFile file("StraightShower.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 1e19);
        Reconstructor::Result result = monte_carlo->RunSingleShower(shower, "straight_shower");
        cout << "Energy, " << shower.Header() << ", " << result.Header() << endl;
        cout << shower.EnergyeV() << ", " << shower.ToString(FriendGroundPlane()) << ", "
             << result.ToString(FriendGroundPlane()) << endl;
    }

    /*
     * Simulate and attempt reconstruction on a typical shower at 10km.
     */
    TEST_F(SampleEvents, TypicalShower)
    {
        TFile file("TypicalShower.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 1, -3), 1e6, -0.1, 1e19);
        Reconstructor::Result result = monte_carlo->RunSingleShower(shower, "typical_shower");
        cout << "Energy, " << shower.Header() << ", " << result.Header() << endl;
        cout << shower.EnergyeV() << ", " << shower.ToString(FriendGroundPlane()) << ", "
             << result.ToString(FriendGroundPlane()) << endl;
    }

    /*
     * Simulate and attempt reconstruction on a straight vertical shower at 30km.
     */
    TEST_F(SampleEvents, DistantShower)
    {
        TFile file("DistantShower.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 3e6, 0, 1e19);
        Reconstructor::Result result = monte_carlo->RunSingleShower(shower, "distant_shower");
        cout <<  endl << "Energy," << shower.Header() << "," << result.Header() << endl;
        cout << shower.EnergyeV() << "," << shower.ToString(FriendGroundPlane()) << ","
             << result.ToString(FriendGroundPlane()) << endl;
    }
}