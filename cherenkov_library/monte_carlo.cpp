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
    void MonteCarlo::ParseFile(boost::property_tree::ptree config)
    {
        // The distribution of shower energies
        std::string energy_formula = "x^(-" + std::to_string(config.get<double>("energy_pow")) + ")";
        energy_distribution = TF1("energy", energy_formula.c_str(), config.get<double>("e_min"),
                                  config.get<double>("e_max"));

        // The distribution of shower vertical directions
        cosine_distribution = TF1("cosine", "cos(x)", -TMath::Pi() / 2, TMath::Pi() / 2);

        // The Monte Carlo distribution of impact parameters
        impact_distribution = TF1("impact", "x", config.get<double>("impact_min"),
                                  config.get<double>("impact_max"));

        // The distribution of first interaction depths (See AbuZayyad 6.1)
        std::string interaction_formula = "e^(-x/ " + std::to_string(config.get<double>("avg_interact")) + ")";
        interact_distribution = TF1("interact", interaction_formula.c_str(), 0, TMath::Infinity());

        // Parameters defining properties of the atmosphere
        scale_height = config.get<double>("atmosphere.scale_height");
        double rho_sea = config.get<double>("atmosphere.rho_sea");
        double detect_elevation = config.get<double>("atmosphere.detect_elevation");
        rho_0 = rho_sea * Exp(-detect_elevation / scale_height);

        // 1 - the index of refraction at the detector (proportional to atmospheric density)
        double delta_sea = 1.0 - config.get<double>("optics.n_air");
        delta_0 = delta_sea * Exp(-detect_elevation / scale_height);
    }

    Shower MonteCarlo::GenerateRandomShower()
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
        TVector3 shower_axis = TVector3(cos(theta), sin(theta) * cos(phi_shower), sin(theta) * sin(phi_shower));

        // We define the origin of both the world and detector frames to be the detector's center of curvature for
        // simplicity. We know that, at the impact point, the position vector of the shower is normal to its direction
        // vector. Let's find a random vector normal to the shower axis.
        TVector3 impact_point = impact_param * RandomPerpendicularVector(shower_axis, rng);

        // Find the depth of the first interaction, the depth of the maximum, and the size of the shower maximum (See
        // AbuZayyad 6.1-6.4). We assume a proton primary.
        // TODO: Add specific numerical values to the configuration file.
        double x_0 = interact_distribution.GetRandom();
        double x_max = 725.0 + 55.0 * (Log(energy) - 18.0) - 70 + x_0;
        double n_max = energy / n_max_ratio;

        // Trace the path of the shower back to the location of the first interaction. Start by finding the elevation of
        // the first interaction.
        double interaction_height = -scale_height * cos(theta) * Log(x_0 / (rho_0 * scale_height * cos(theta)));
        double param = (interaction_height - impact_point.Z()) / (shower_axis.Z());
        TVector3 starting_position = impact_point + param * shower_axis;

        // Create a new shower with all of the randomly determined parameters.
        Shower::Params params;
        params.x_0 = x_0;
        params.x_max = x_max;
        params.n_max = n_max;
        params.rho_0 = rho_0;
        params.scale_height = scale_height;
        params.delta_0 = delta_0;
        return Shower(params, starting_position, shower_axis);
    }
    
}
