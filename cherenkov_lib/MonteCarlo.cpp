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

using namespace std;
using namespace boost::property_tree;
using namespace TMath;

namespace cherenkov_simulator
{
    MonteCarlo::MonteCarlo(const ptree& config) : simulator(config), reconstructor(config)
    {
        n_showers = config.get<int>("simulation.n_showers");
        elevation = config.get<double>("surroundings.elevation");

        energy_pow = config.get<double>("monte_carlo.energy_pow");
        energy_min = config.get<double>("monte_carlo.energy_min");
        energy_max = config.get<double>("monte_carlo.energy_max");
        impact_min = config.get<double>("monte_carlo.impact_min");
        impact_max = config.get<double>("monte_carlo.impact_max");
        begn_depth = config.get<double>("monte_carlo.begn_depth");
    }

    void MonteCarlo::PerformMonteCarlo(string output_file) const
    {
        TFile file((output_file + ".root").c_str(), "RECREATE");
        ofstream fout = ofstream(output_file + ".csv");
        unsigned int start_seed = gRandom->GetSeed();
        fout << "Seed, ID, Energy, " << Shower::Header() << ", " << Reconstructor::Result::Header() << endl;

        Plane ground_plane = simulator.GroundPlane();
        for (int i = 0; i < n_showers;)
        {
            Shower shower = GenerateShower();
            Reconstructor::Result result = RunSingleShower(shower, to_string(i));
            if (!result.triggered) continue;
            else i++;
            cout << "Shower " << i << " finished" << endl;
            fout << start_seed << ", " << i << ", " << shower.EnergyeV() << ", " << shower.ToString(ground_plane)
                 << ", " << result.ToString(ground_plane) << endl;
        }
    }

    Reconstructor::Result MonteCarlo::RunSingleShower(Shower shower, string ident) const
    {
        PhotonCount data;
        try
        {
            data = simulator.SimulateShower(shower);
        }
        catch (out_of_range& err)
        {
            cout << err.what() << endl;
            cout << "Skipping this shower..." << endl;
            return Reconstructor::Result();
        }
        if (data.Empty())
        {
            Reconstructor::Result result = Reconstructor::Result();
            result.triggered = false;
            return result;
        }

        TH2I befor_noise_pixl = Analysis::MakePixlProfile(data, ident + "_befor_noise_pixl");
        TGraph befor_noise_time = Analysis::MakeTimeProfile(data);
        reconstructor.AddNoise(data);
        TH2I after_noise_pixl = Analysis::MakePixlProfile(data, ident + "_after_noise_pixl");
        TGraph after_noise_time = Analysis::MakeTimeProfile(data);
        reconstructor.ClearNoise(data);
        TH2I after_clear_pixl = Analysis::MakePixlProfile(data, ident + "_after_clear_pixl");
        TGraph after_clear_time = Analysis::MakeTimeProfile(data);

        Reconstructor::Result result = reconstructor.Reconstruct(data);
        if (!result.triggered) return result;

        befor_noise_pixl.Write();
        befor_noise_time.Write((ident + "_befor_noise_time").c_str());
        after_noise_pixl.Write();
        after_noise_time.Write((ident + "_after_noise_time").c_str());
        after_clear_pixl.Write();
        after_clear_time.Write((ident + "_after_clear_time").c_str());

        Plane ground_plane = simulator.GroundPlane();
        shower.Direction().Write((ident + "_orig_direction").c_str());
        shower.PlaneImpact(ground_plane).Write((ident + "_orig_gnd_impact").c_str());
        result.mono_recon.Direction().Write((ident + "_mono_direction").c_str());
        result.mono_recon.PlaneImpact(ground_plane).Write((ident + "_mono_gnd_impact").c_str());
        result.chkv_recon.Direction().Write((ident + "_chkv_direction").c_str());
        result.chkv_recon.PlaneImpact(ground_plane).Write((ident + "_chkv_gnd_impact").c_str());
        return result;
    }

    Shower MonteCarlo::GenerateShower() const
    {
        double zenith = Utility::RandCosine();
        double azmuth = gRandom->Uniform(TwoPi());
        TVector3 axis = TVector3(sin(zenith) * cos(azmuth), sin(zenith) * sin(azmuth), -cos(zenith));

        double im_par = Utility::RandLinear(impact_min, impact_max);
        double im_ang = gRandom->Uniform(TwoPi());
        double energy = Utility::RandPower(energy_min, energy_max, energy_pow);
        return GenerateShower(axis, im_par, im_ang, energy);
    }

    Shower MonteCarlo::GenerateShower(TVector3 axis, double im_par, double im_ang, double energy) const
    {
        // Start with an impact point directly in front of the detector, then rotate it by a random angle.
        TVector3 impact_pos = TVector3(1, 0, 0).Cross(axis).Unit();
        impact_pos.Rotate(im_ang, axis);
        impact_pos *= im_par;

        double start_h = -scale_h * Log(begn_depth * Abs(axis.CosTheta()) / (rho_sea * scale_h)) - elevation;
        double trace = (start_h - impact_pos.Z()) / (axis.Z());
        TVector3 start_pos = impact_pos + trace * axis;
        return Shower(energy, elevation, start_pos, axis);
    }

    int MonteCarlo::Run(int argc, const char* argv[])
    {
        string output_file = "Output";
        string config_file = "Config.xml";
        if (argc > 1) output_file = string(argv[1]);
        if (argc > 2) config_file = string(argv[2]);
        try
        {
            ptree config = Utility::ParseXMLFile(config_file).get_child("config");
            if (config.get<bool>("simulation.time_seed")) gRandom->SetSeed();
            if (argc > 3) gRandom->SetSeed(stoul(argv[3]));
            MonteCarlo(config).PerformMonteCarlo(output_file);
            return 0;
        }
        catch (runtime_error& err)
        {
            cout << err.what() << endl;
            return -1;
        }
    }
}
