// MonteCarlo.cpp
//
// Author: Matthew Dutson
//
// Implementation of MonteCarlo.h

#include <TMath.h>

#include "MonteCarlo.h"
#include "Utility.h"

using namespace TMath;

using std::cout;
using std::endl;

namespace cherenkov_simulator
{
    MonteCarlo::MonteCarlo(boost::property_tree::ptree config): simulator(config), reconstructor(config)
    {
        // The distribution of shower energies
        std::string energy_formula = "x^(-" + std::to_string(energy_pow) + ")";
        double e_min = config.get<double>("monte_carlo.e_min");
        double e_max = config.get<double>("monte_carlo.e_max");
        energy_distribution = TF1("energy", energy_formula.c_str(), e_min, e_max);

        // The distribution of shower vertical directions
        cosine_distribution = TF1("cosine", "cos(x)", 0.0, TMath::Pi() / 2);

        // The Monte Carlo distribution of impact parameters
        double impact_min = config.get<double>("monte_carlo.impact_min");
        double impact_max = config.get<double>("monte_carlo.impact_max");
        impact_distribution = TF1("impact", "x", impact_min, impact_max);

        // First interaction depths follow an exponential distribution (See AbuZayyad 6.1)
        start_tracking = config.get<double>("monte_carlo.start_tracking");

        // Determine the local atmosphere based on the elevation
        double detect_elevation = config.get<double>("surroundings.detect_elevation");
        rho_0 = rho_sea * Exp(-detect_elevation / scale_height);
        delta_0 = (refrac_sea - 1) * Exp(-detect_elevation / scale_height);

        // The number of showers to simulate in PerformMonteCarlo
        n_showers = config.get<int>("simulation.n_showers");

        // The random number generator
        rng = TRandom3();
    }

    void MonteCarlo::PerformMonteCarlo()
    {
        cout << "Running simulation";
        int n_triggered = 0;
        int n_ground_used = 0;
        double mono_rp_err = 0.0, ckv_rp_err = 0.0;
        double mono_psi_err = 0.0, ckv_psi_err = 0.0;
        for (int i = 0; i < n_showers; i++)
        {
            cout << ".";
            Shower shower = GenerateShower();
            PhotonCount data = simulator.SimulateShower(shower);
            reconstructor.AddNoise(&data);
            bool triggered, ground_used;
            Shower mono_shower = reconstructor.Reconstruct(data, false, &triggered, &ground_used);
            Shower ckv_shower = reconstructor.Reconstruct(data, true, &triggered, &ground_used);
            if (triggered)
            {
                n_triggered++;
                mono_rp_err += Utility::PercentError(shower.Position().Mag(), mono_shower.Position().Mag());
                mono_psi_err += shower.Direction().Angle(mono_shower.Direction());
            }
            if (triggered && ground_used)
            {
                n_ground_used++;
                ckv_rp_err += Utility::PercentError(shower.Position().Mag(), ckv_shower.Position().Mag());
                ckv_psi_err += shower.Direction().Angle(ckv_shower.Direction());
            }
        }
        mono_psi_err /= (double) n_triggered;
        ckv_psi_err /= (double) n_ground_used;
        mono_rp_err /= (double) n_triggered;
        ckv_rp_err /= (double) n_ground_used;

        cout << endl << "Number simulated: \t" << n_showers << endl;
        cout << "Number triggered: \t" << n_triggered << endl;
        cout << "Number with ground point found: \t" << n_ground_used << endl;
        cout << "Mean monocular impact error: \t" << mono_rp_err << " %" << endl;
        cout << "Mean Cherenkov impact error: \t" << ckv_rp_err << " %" << endl;
        cout << "Mean monocular angle error: \t" << mono_psi_err << " rad" << endl;
        cout << "Mean Cherenkov angle error: \t" << ckv_psi_err << " rad" << endl;
    }

    Shower MonteCarlo::GenerateShower()
    {
        // Determine the direction of the shower and its position relative to the detector. The angle of the shower
        // relative to the vertical goes as cos(theta) because shower have an isotropic flux in space.
        double theta = cosine_distribution.GetRandom();
        double phi_shower = rng.Uniform(TwoPi());

        // Determine the impact parameter.
        double impact_param = impact_distribution.GetRandom();

        // Find the Cartesian shower axis vector. This vector is in the world frame (z is normal to the surface of the
        // earth, with x and y parallel to the surface. Note that the surface of the earth may not be parallel to the
        // local ground.
        TVector3 shower_axis = TVector3(sin(theta) * cos(phi_shower), sin(theta) * sin(phi_shower), -cos(theta));
        return GenerateShower(shower_axis, impact_param, rng.Uniform(TwoPi()), energy_distribution.GetRandom());
    }

    Shower MonteCarlo::GenerateShower(TVector3 axis, double impact_param, double impact_angle, double energy)
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
        Shower::Params params;
        params.energy = energy;
        params.x_max = x_max;
        params.n_max = n_max;
        params.rho_0 = rho_0;
        params.scale_height = scale_height;
        params.delta_0 = delta_0;
        return Shower(params, starting_position, axis);
    }
}
