//
// Created by Matthew Dutson on 2/7/17.
//

#include <gtest/gtest.h>
#include <TFile.h>
#include <TVectorD.h>
#include "reconstructor.h"
#include "simulator.h"
#include "monte_carlo.h"
#include "utility.h"

using boost::property_tree::ptree;

using namespace cherenkov_library;

namespace cherenkov_tests
{
    TEST(reconstructor, straight_down)
    {
        // Parse the XML file and give it to the simulator and monte carlo.
        ptree config = ParseXMLFile("../../config.xml").get_child("config");
        Simulator simulator = Simulator(config);
        MonteCarlo monte_carlo = MonteCarlo(config);
        Reconstructor reconstructor = Reconstructor(config);

        // Construct and simulate a shower at 10 km. The axis coordinates are in the world frame.
        Shower shower = monte_carlo.GenerateShower(TVector3(0, 0, -1), 1000000, 0, 10e19);
        PhotonCount data = simulator.SimulateShower(shower);

        // Attempt to reconstruct the shower plane and geometry.
        Plane plane = reconstructor.FitSDPlane(data);
        double t_0, impact, angle;
        reconstructor.TimeProfileFit(data, plane, &t_0, &impact, &angle);

        // Write results to a file.
        TFile file("../../reconstructor_straight_down.root", "RECREATE");
        TVectorD reconstruction(3);
        reconstruction[0] = t_0;
        reconstruction[1] = impact;
        reconstruction[2] = angle;
        reconstruction.Write("reconstruction");
    }
}