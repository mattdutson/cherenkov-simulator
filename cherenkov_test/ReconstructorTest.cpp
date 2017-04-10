// ReconstructorTest.cpp
//
// Author: Matthew Dutson
//
// Tests of Reconstructr.h

#include <gtest/gtest.h>
#include <TFile.h>

#include "Reconstructor.h"
#include "Simulator.h"
#include "MonteCarlo.h"
#include "Utility.h"
#include "Analysis.h"

using namespace TMath;

using boost::property_tree::ptree;
using std::vector;
using std::string;

namespace cherenkov_simulator
{
    class ReconstructorTest : public testing::Test
    {

    protected:

        Simulator* simulator;
        MonteCarlo* monte_carlo;
        Reconstructor* reconstructor;

        virtual void SetUp()
        {
            ptree config = Utility::ParseXMLFile("../../../Config.xml").get_child("config");
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

        Bool1D FriendTriggeringState(PhotonCount& data)
        {
            return reconstructor->GetTriggeringState(data);
        }

        Bool3D FriendSixSigma(PhotonCount& data)
        {
            return reconstructor->GetThresholdMatrices(data, 6, true);
        }

        TRotation FriendFitSDPlane(PhotonCount& data)
        {
            return reconstructor->FitSDPlane(data);
        }

        Shower FriendMonocularFit(PhotonCount& data, TRotation to_sd_plane, string graph_file = "")
        {
            return reconstructor->MonocularFit(data, to_sd_plane, graph_file);
        }

        bool FriendFindGroundImpact(PhotonCount& data, TVector3& impact)
        {
            return reconstructor->FindGroundImpact(data, impact);
        }
    };

    TEST_F(ReconstructorTest, StraightShower)
    {
        // Construct and simulate a shower at 10 km. The axis coordinates are in the world frame
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);
        reconstructor->AddNoise(data);
        reconstructor->ClearNoise(data);

        // Attempt to reconstruct the shower plane and geometry, writing the fit graph to a file
        TRotation to_sd_plane = FriendFitSDPlane(data);
        FriendMonocularFit(data, to_sd_plane, "StraightShowerRecon.root");
    }

    TEST_F(ReconstructorTest, AngleShower)
    {
        // Construct and simulate a shower at 10km which is skewed at some angle.
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 0, -2), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);
        reconstructor->AddNoise(data);
        reconstructor->ClearNoise(data);

        // Attempt to reconstruct the shower plane and geometry, writing the fit graph to a file
        TRotation to_sd_plane = FriendFitSDPlane(data);
        FriendMonocularFit(data, to_sd_plane, "AngleShowerRecon.root");
    }

    TEST_F(ReconstructorTest, AddSubtractNoise)
    {
        TFile file("AddSubtractNoise.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 1, -3), 1e6, -0.1, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        Analysis::MakeProfileGraph(data).Write("before_noise_graph");
        Analysis::MakeSumMap(data).Write("before_noise_map");

        reconstructor->AddNoise(data);
        Analysis::MakeProfileGraph(data).Write("after_noise_graph");
        Analysis::MakeSumMap(data).Write("after_noise_map");

        reconstructor->ClearNoise(data);
        Analysis::MakeProfileGraph(data).Write("after_clear_graph");
        Analysis::MakeSumMap(data).Write("after_clear_map");
    }

    TEST_F(ReconstructorTest, TriggeringMaps)
    {
        TFile file("TriggeringMaps.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        Bool3D triggering_matrices = FriendSixSigma(data);
        Bool1D triggering_state = FriendTriggeringState(data);
        for (int i = 0; i < triggering_state.size(); i++)
        {
            TH2C frame_map = Analysis::GetBooleanMap(triggering_matrices[i]);
            std::string write_name;
            if (triggering_state[i])
            {
                write_name = std::to_string(i) + "tr";
            }
            else
            {
                write_name = std::to_string(i) + "nt";
            }
            frame_map.Write(write_name.c_str());
        }
    }

    TEST_F(ReconstructorTest, FindGroundImpact)
    {
        TFile file("ImpactPoint.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        TVector3 impact;
        if (FriendFindGroundImpact(data, impact))
        {
            impact.Write("ground_impact");
        }
    }

    TEST_F(ReconstructorTest, StraightShowerGeometry)
    {
        TFile file("StraightGeometry.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        Reconstructor::Result result = reconstructor->Reconstruct(data);
        result.mono.Position().Write("straight_shower_position");
        result.mono.Direction().Write("straight_shower_direction");
    }

    TEST_F(ReconstructorTest, AngleShowerGeometry)
    {
        TFile file("AngleGeometry.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 0, -2), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        Reconstructor::Result result = reconstructor->Reconstruct(data);
        result.mono.Position().Write("angle_shower_position");
        result.mono.Direction().Write("angle_shower_direction");
    }

    TEST_F(ReconstructorTest, StraightCherenkovRecon)
    {
        TFile file("StraightCherenkovRecon.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);
        reconstructor->AddNoise(data);

        Reconstructor::Result result = reconstructor->Reconstruct(data);
        result.ckv.Position().Write("straight_shower_position");
        result.ckv.Direction().Write("straight_shower_direction");
        std::cout << "Mono Direction Error: " << result.mono.Direction().Angle(shower.Direction()) << std::endl;
        std::cout << "Mono Impact Error: " << Abs(result.mono.Position().Mag() - 1e6) << std::endl;
        std::cout << "Cherenkov Direction Error: " << result.ckv.Direction().Angle(shower.Direction()) << std::endl;
        std::cout << "Cherenkov Impact Error: " << Abs(result.ckv.Position().Mag() - 1e6) << std::endl;
    }

    TEST_F(ReconstructorTest, AngleCherenkovRecon)
    {
        TFile file("AngleCherenkovRecon.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 0, -2), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);
        reconstructor->AddNoise(data);

        Reconstructor::Result result = reconstructor->Reconstruct(data);
        result.ckv.Position().Write("angle_shower_position");
        result.ckv.Direction().Write("angle_shower_direction");
        std::cout << "Mono Direction Error: " << result.mono.Direction().Angle(shower.Direction()) << std::endl;
        std::cout << "Mono Impact Error: " << Abs(result.mono.Position().Mag() - 1e6) << std::endl;
        std::cout << "Cherenkov Direction Error: " << result.ckv.Direction().Angle(shower.Direction()) << std::endl;
        std::cout << "Cherenkov Impact Error: " << Abs(result.ckv.Position().Mag() - 1e6) << std::endl;
    }

    TEST_F(ReconstructorTest, StackOverflow)
    {
        // Rebuild shower 97 (i = 1196)
        Shower::Params params = Shower::Params();
        params.energy = 9.7224104168849578e20;
        params.x_max= 889.32756726164371;
        params.n_max = 699453986826.25598;
        params.rho_0 = 0.0010354826266577842;
        params.scale_height = 841300;
        params.delta_0 = 0.0002451346626372697;
        TVector3 position = TVector3(-405043.23971585685, -3832795.7567525646, 3526742.1184894536);
        TVector3 direction = TVector3(471764977.74197131, 28917762299.766682, -7895622925.3915901);
        Shower shower = Shower(params, position, direction);

        // Simulate and attempt reconstruction on the shower. There may be a memory access error (caused by a stack
        // overflow) during the ClearNoise step.
        PhotonCount data = simulator->SimulateShower(shower);
        reconstructor->AddNoise(data);
        reconstructor->ClearNoise(data);
        reconstructor->Reconstruct(data);
    }
}