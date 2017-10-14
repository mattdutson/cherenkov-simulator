// MonteCarlo.cpp
//
// Author: Matthew Dutson
//
// Implementation of MonteCarlo.h

#include <fstream>
#include <TMath.h>
#include <TFile.h>

#include "MonteCarlo.h"
#include "Analysis.h"

using namespace TMath;

using std::cout;
using std::endl;
using boost::property_tree::ptree;

namespace cherenkov_simulator
{
    MonteCarlo::MonteCarlo(const boost::property_tree::ptree& config) : simulator(config), reconstructor(config)
    {
        // The distribution of shower energies
        energy_pow = config.get<double>("monte_carlo.energy_pow");
        e_min = config.get<double>("monte_carlo.e_min");
        e_max = config.get<double>("monte_carlo.e_max");

        // The Monte Carlo distribution of impact parameters
        impact_min = config.get<double>("monte_carlo.impact_min");
        impact_max = config.get<double>("monte_carlo.impact_max");

        // First interaction depths follow an exponential distribution (See AbuZayyad 6.1)
        start_tracking = config.get<double>("monte_carlo.start_tracking");

        // Determine the local atmosphere based on the elevation
        auto detect_elevation = config.get<double>("surroundings.detect_elevation");
        rho_0 = rho_sea * Exp(-detect_elevation / scale_height);
        delta_0 = (refrac_sea - 1) * Exp(-detect_elevation / scale_height);

        // The number of showers to simulate in PerformMonteCarlo
        n_showers = config.get<int>("simulation.n_showers");
    }

    void MonteCarlo::PerformMonteCarlo(std::string out_file) const
    {
        // Open a ROOT file to store information about each shower
        TFile file((out_file + ".root").c_str(), "RECREATE");
        std::ofstream fout = std::ofstream(out_file + ".csv");
        unsigned int start_seed = gRandom->GetSeed();
        fout << "Seed, ID, Energy, " << Shower::Header() << ", " << Reconstructor::Result::Header() << endl;

        // Simulate a user-defined number of showers
        for (int i = 0; i < n_showers;)
        {
            // Simulate the shower and record photon counts before noise is added
            Shower shower = GenerateShower();
            PhotonCount data;
            try
            {
                data = simulator.SimulateShower(shower);
            }
            catch (std::out_of_range& err)
            {
                cout << err.what() << endl;
                cout << "Skipping this shower..." << endl;
                continue;
            }
            if (data.Empty()) continue;
            else i++;

            // Write a map of the initial shower track
            Analysis::MakeSumMap(data).Write((std::to_string(i) + "_before_noise_map").c_str());
            Analysis::MakeProfileGraph(data).Write((std::to_string(i) + "_before_noise_graph").c_str());

            // Add noise and record the new photon counts
            reconstructor.AddNoise(data);
            Analysis::MakeSumMap(data).Write((std::to_string(i) + "_after_noise_map").c_str());
            Analysis::MakeProfileGraph(data).Write((std::to_string(i) + "_after_noise_graph").c_str());

            // Clear noise and record the new photon counts
            reconstructor.ClearNoise(data);
            Analysis::MakeSumMap(data).Write((std::to_string(i) + "_after_clear_map").c_str());
            Analysis::MakeProfileGraph(data).Write((std::to_string(i) + "_after_clear_graph").c_str());

            // Attempt both monocular and hybrid reconstruction of the shower
            Reconstructor::Result result = reconstructor.Reconstruct(data);
            cout << "Shower " << i << " finished" << endl;
            fout << start_seed << ", " << i << ", " << shower.EnergyeV() << ", " << shower.ToString() << ", "
                 << result.ToString() << endl;
        }
    }

    Shower MonteCarlo::GenerateShower() const
    {
        // Determine the direction of the shower and its position relative to the detector. The angle of the shower
        // relative to the vertical goes as cos(theta) because shower have an isotropic flux in space.
        double theta = Utility::RandCosine();
        double phi_shower = gRandom->Uniform(TwoPi());
        double energy = Utility::RandPower(e_min, e_max, energy_pow);

        // Determine the impact parameter.
        double impact_param = Utility::RandLinear(impact_min, impact_max);

        // Find the Cartesian shower axis vector. This vector is in the world frame (z is normal to the surface of the
        // earth, with x and y parallel to the surface. Note that the surface of the earth may not be parallel to the
        // local ground.
        TVector3 shower_axis = TVector3(sin(theta) * cos(phi_shower), sin(theta) * sin(phi_shower), -cos(theta));
        return GenerateShower(shower_axis, impact_param, gRandom->Uniform(TwoPi()), energy);
    }

    Shower MonteCarlo::GenerateShower(TVector3 axis, double impact_param, double impact_angle, double energy) const
    {
        // We define the origin of both the world and detector frames to be the detector's center of curvature for
        // simplicity. We know that, at the impact point, the position vector of the shower is normal to its direction
        // vector.
        TVector3 impact_direction = TVector3(1, 0, 0).Cross(axis).Unit();
        impact_direction.Rotate(impact_angle, axis);
        TVector3 impact_point = impact_param * impact_direction;

        // Find the depth of the first interaction, the depth of the maximum, and the size of the shower maximum (See
        // AbuZayyad 6.1-6.4). We assume a proton primary.
        double x_max = x_max_1 + x_max_2 * (Log10(energy) - x_max_3);
        double n_max = energy / n_max_ratio;

        // Trace the path of the shower back to the location of the first interaction. Start by finding the elevation of
        // the first interaction. See notes from 1/25.
        double cos_theta = Abs(axis.CosTheta());
        double interaction_height = -scale_height * Log(start_tracking * cos_theta / (rho_0 * scale_height));
        double param = (interaction_height - impact_point.Z()) / (axis.Z());
        TVector3 starting_position = impact_point + param * axis;

        // Create a new shower with all of the randomly determined parameters.
        Shower::Params params = Shower::Params();
        params.energy = energy;
        params.x_max = x_max;
        params.n_max = n_max;
        params.rho_0 = rho_0;
        params.scale_height = scale_height;
        params.delta_0 = delta_0;
        return Shower(params, starting_position, axis);
    }

    int MonteCarlo::Run(int argc, const char* argv[])
    {
        std::string out_file = "Output";
        std::string config_file = "Config.xml";
        if (argc > 1) out_file = std::string(argv[1]);
        if (argc > 2) config_file = std::string(argv[2]);
        try
        {
            ptree config = Utility::ParseXMLFile(config_file).get_child("config");
            if (config.get<bool>("simulation.time_seed")) gRandom->SetSeed();
            if (argc > 3) gRandom->SetSeed(std::stoul(argv[3]));
            MonteCarlo(config).PerformMonteCarlo(out_file);
            return 0;
        }
        catch (std::runtime_error& err)
        {
            cout << err.what() << endl;
            return -1;
        }
    }
}
