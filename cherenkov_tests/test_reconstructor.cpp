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
#include "data_analysis.h"

using boost::property_tree::ptree;

using namespace cherenkov_library;

namespace cherenkov_tests
{
    class ReconstructorTest : public testing::Test
    {

    protected:

        Simulator* simulator;
        MonteCarlo* monte_carlo;
        Reconstructor* reconstructor;

        virtual void SetUp()
        {
            ptree config = ParseXMLFile("../../config.xml").get_child("config");
            simulator = new Simulator(config);
            monte_carlo = new MonteCarlo(config);
            reconstructor = new Reconstructor(config);
        }

        virtual void TearDown()
        {
            delete simulator;
            delete monte_carlo;
            delete reconstructor;
        }
    };

    TEST_F(ReconstructorTest, StraightShower)
    {
        // Construct and simulate a shower at 10 km. The axis coordinates are in the world frame.
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 10e19);
        PhotonCount data = simulator->SimulateShower(shower);

        // Attempt to reconstruct the shower plane and geometry.
        Plane plane = reconstructor->FitSDPlane(data);
        double t_0, impact, angle;
        TGraph time_profile = reconstructor->MonocularFit(data, plane, &t_0, &impact, &angle);

        // Write results to a file.
        TFile file("../../cherenkov_tests/straight_shower_recon.root", "RECREATE");
        time_profile.Write("straight_shower_profile");
        TVectorD reconstruction(3);
        reconstruction[0] = t_0;
        reconstruction[1] = impact;
        reconstruction[2] = angle;
        reconstruction.Write("straight_shower_params");
    }

    TEST_F(ReconstructorTest, AngleShower)
    {
        // Construct and simulate a shower at 10km which is skewed at some angle.
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 0, -2), 1e6, 0, 10e19);
        PhotonCount data = simulator->SimulateShower(shower);

        // Attempt to reconstruct the shower plane and geometry.
        Plane plane = reconstructor->FitSDPlane(data);
        double t_0, impact, angle;
        TGraph time_profile = reconstructor->MonocularFit(data, plane, &t_0, &impact, &angle);

        // Draw a map of impacts.
        TFile file("../../cherenkov_tests/angle_shower_recon.root", "RECREATE");
        time_profile.Write("angle_shower_profile");
        TVectorD reconstruction(3);
        reconstruction[0] = t_0;
        reconstruction[1] = impact;
        reconstruction[2] = angle;
        reconstruction.Write("angle_shower_params");
    }

    TEST_F(ReconstructorTest, AddSubtractNoise)
    {
        TFile file("../../cherenkov_tests/add_subtract_noise.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 10e19);
        PhotonCount data = simulator->SimulateShower(shower);

        DataAnalysis::MakeProfileGraph(data).Write("before_noise_graph");
        DataAnalysis::MakeSumMap(data).Write("before_noise_map");

        simulator->AddNoise(&data);
        DataAnalysis::MakeProfileGraph(data).Write("after_noise_graph");
        DataAnalysis::MakeSumMap(data).Write("after_noise_map");

        reconstructor->SubtractAverageNoise(&data);
        DataAnalysis::MakeProfileGraph(data).Write("after_subtract_graph");
        DataAnalysis::MakeSumMap(data).Write("after_subtract_map");

        reconstructor->ApplyTriggering(&data);
        DataAnalysis::MakeProfileGraph(data).Write("after_trigger_graph");
        DataAnalysis::MakeSumMap(data).Write("after_trigger_map");
    }

    TEST_F(ReconstructorTest, TriggeringMaps)
    {
        TFile file("../../cherenkov_tests/triggering_maps.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 10e19);
        PhotonCount data = simulator->SimulateShower(shower);

        std::vector<std::vector<std::vector<bool>>> triggering_matrices = reconstructor->GetTriggeringMatrices(data);
        for (int i = 0; i < triggering_matrices.size(); i++)
        {
            TH2C frame_map = DataAnalysis::GetBooleanMap(triggering_matrices[i]);
            std::string write_name;
            if (reconstructor->FrameTriggered(triggering_matrices[i]))
            {
                write_name = std::to_string(i) + "_tr";
            }
            else
            {
                write_name = std::to_string(i) + "_nt";
            }
            frame_map.Write(write_name.c_str());
        }
    }

    TEST_F(ReconstructorTest, FindGroundImpact)
    {
        TFile file("../../cherenkov_tests/impact_point.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 10e19);
        PhotonCount data = simulator->SimulateShower(shower);

        TVector3 impact;
        if (reconstructor->FindGroundImpact(data, &impact))
        {
            impact.Write("ground_impact");
        }
    }
}