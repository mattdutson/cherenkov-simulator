// ReconstructorTest.cpp
//
// Author: Matthew Dutson
//
// Tests of Reconstructr.h

#include <gtest/gtest.h>
#include <TFile.h>
#include <TVectorD.h>

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

        typedef std::vector<std::vector<std::vector<bool>>> Bool3D;

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

        bool FriendFrameTriggered(int i, Bool3D* triggers) {return reconstructor->FrameTriggered(i, triggers);}

        Bool3D FriendTriggeringMatrices(PhotonCount data) {return reconstructor->GetThresholdMatrices(data, 0);}

        TRotation FriendFitSDPlane(PhotonCount data) {return reconstructor->FitSDPlane(data);}

        Shower FriendMonocularFit(PhotonCount data, TRotation to_sd_plane, string graph_file = "") {return reconstructor->MonocularFit(data, to_sd_plane, graph_file);}

        Shower FriendHybridFit(PhotonCount data, TVector3 impact, TRotation to_sd_plane, string graph_file = "") {return reconstructor->HybridFit(data, impact, to_sd_plane, graph_file);}

        void FriendSubtractAverageNoise(PhotonCount* data) {reconstructor->SubtractAverageNoise(data);}

        void FriendThreeSigmaFilter(PhotonCount* data) {reconstructor->ThreeSigmaFilter(data);}

        void FriendRecursiveSearch(PhotonCount* data) {reconstructor->RecursiveSearch(data);}

        bool FriendFindGroundImpact(PhotonCount data, TVector3* impact) {return reconstructor->FindGroundImpact(data, impact);}
    };

    TEST_F(ReconstructorTest, StraightShower)
    {
        // Construct and simulate a shower at 10 km. The axis coordinates are in the world frame
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);
        reconstructor->AddNoise(&data);
        reconstructor->ClearNoise(&data);

        // Attempt to reconstruct the shower plane and geometry, writing the fit graph to a file
        TRotation to_sd_plane = FriendFitSDPlane(data);
        FriendMonocularFit(data, to_sd_plane, "StraightShowerRecon.root");
    }

    TEST_F(ReconstructorTest, AngleShower)
    {
        // Construct and simulate a shower at 10km which is skewed at some angle.
        Shower shower = monte_carlo->GenerateShower(TVector3(1, 0, -2), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);
        reconstructor->AddNoise(&data);
        reconstructor->ClearNoise(&data);

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

        reconstructor->AddNoise(&data);
        Analysis::MakeProfileGraph(data).Write("after_noise_graph");
        Analysis::MakeSumMap(data).Write("after_noise_map");

        reconstructor->ClearNoise(&data);
        Analysis::MakeProfileGraph(data).Write("after_clear_graph");
        Analysis::MakeSumMap(data).Write("after_clear_map");
    }

    TEST_F(ReconstructorTest, TriggeringMaps)
    {
        TFile file("TriggeringMaps.root", "RECREATE");
        Shower shower = monte_carlo->GenerateShower(TVector3(0, 0, -1), 1e6, 0, 1e19);
        PhotonCount data = simulator->SimulateShower(shower);

        std::vector<std::vector<std::vector<bool>>> triggering_matrices = FriendTriggeringMatrices(data);
        for (int i = 0; i < triggering_matrices.size(); i++)
        {
            TH2C frame_map = Analysis::GetBooleanMap(triggering_matrices[i]);
            std::string write_name;
            if (FriendFrameTriggered(i, &triggering_matrices))
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
        if (FriendFindGroundImpact(data, &impact))
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
        reconstructor->AddNoise(&data);

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
        reconstructor->AddNoise(&data);

        Reconstructor::Result result = reconstructor->Reconstruct(data);
        result.ckv.Position().Write("angle_shower_position");
        result.ckv.Direction().Write("angle_shower_direction");
        std::cout << "Mono Direction Error: " << result.mono.Direction().Angle(shower.Direction()) << std::endl;
        std::cout << "Mono Impact Error: " << Abs(result.mono.Position().Mag() - 1e6) << std::endl;
        std::cout << "Cherenkov Direction Error: " << result.ckv.Direction().Angle(shower.Direction()) << std::endl;
        std::cout << "Cherenkov Impact Error: " << Abs(result.ckv.Position().Mag() - 1e6) << std::endl;
    }
}