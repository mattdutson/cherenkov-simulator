// accuracy_estimator.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include <TMath.h>

#include "monte_carlo.h"
#include "utility.h"

using namespace TMath;

namespace cherenkov_library
{
    MonteCarlo::MonteCarlo(boost::property_tree::ptree config)
    {
        TRotation rotate_to_world = MakeRotation(config.get<double>("elevation_angle"));
        detector_axis = rotate_to_world * TVector3(0, 0, 1);

        // The distribution of shower energies
        std::string energy_formula = "x^(-" + std::to_string(config.get<double>("energy_pow")) + ")";
        energy_distribution = TF1("energy", energy_formula.c_str(), config.get<double>("e_min"),
                                  config.get<double>("e_max"));

        // The distribution of shower vertical directions
        cosine_distribution = TF1("cosine", "cos(x)", -TMath::Pi() / 2, TMath::Pi() / 2);

        // The Monte Carlo distribution of impact parameters
        impact_distribution = TF1("impact", "x", config.get<double>("impact_min"),
                                  config.get<double>("impact_max"));

        // First interaction depths follow an exponential distribution (See AbuZayyad 6.1)
        first_interact = config.get<double>("first_interact");
        start_tracking = config.get<double>("start_tracking");

        // Parameter for determining n_max for a shower.
        n_max_ratio = config.get<double>("n_max_ratio");

        // Parameters defining properties of the atmosphere
        scale_height = config.get<double>("scale_height");
        double rho_sea = config.get<double>("rho_sea");
        double detect_elevation = config.get<double>("detect_elevation");
        rho_0 = rho_sea * Exp(-detect_elevation / scale_height);

        // 1 - the index of refraction at the detector (proportional to atmospheric density)
        double delta_sea = config.get<double>("refrac_sea") - 1.0;
        delta_0 = delta_sea * Exp(-detect_elevation / scale_height);

        // Parameters used when determing the depth of the shower maximum
        x_max_1 = config.get<double>("x_max_1");
        x_max_2 = config.get<double>("x_max_2");
        x_max_3 = config.get<double>("x_max_3");
    }

    Shower MonteCarlo::GenerateShower()
    {
        double energy = energy_distribution.GetRandom();

        // Determine the direction of the shower and its position relative to the detector. The angle of the shower
        // relative to the vertical goes as cos(theta) because shower have an isotropic flux in space.
        double theta = cosine_distribution.GetRandom();
        double phi_shower = rng.Uniform(TwoPi());

        // Determine the impact parameter.
        double impact_param = impact_distribution.GetRandom();

        // Find the Cartesian shower axis vector. This vector is in the world frame (z is normal to the surface of the
        // earth, with x and y parallel to the surface. Note that the surface of the earth may not be parallel to the
        // local ground.
        TVector3 shower_axis = TVector3(sin(theta) * cos(phi_shower), sin(theta) * sin(phi_shower), cos(theta));
        return GenerateShower(shower_axis, impact_param, rng.Uniform(TwoPi()));
    }

    Shower MonteCarlo::GenerateShower(TVector3 axis, double impact_param, double impact_angle)
    {
        return GenerateShower(axis, impact_param, impact_angle, energy_distribution.GetRandom());
    }

    Shower MonteCarlo::GenerateShower(TVector3 axis, double impact_param, double impact_angle, double energy)
    {
        // We define the origin of both the world and detector frames to be the detector's center of curvature for
        // simplicity. We know that, at the impact point, the position vector of the shower is normal to its direction
        // vector.
        TVector3 impact_direction = TVector3(1, 0, 0).Cross(axis);
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
        params.x_0 = first_interact;
        params.x_max = x_max;
        params.n_max = n_max;
        params.rho_0 = rho_0;
        params.scale_height = scale_height;
        params.delta_0 = delta_0;
        return Shower(params, starting_position, axis);
    }
}
