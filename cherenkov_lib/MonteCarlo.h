// monte_carlo.h
// cherenkov_lib
//
// Created by Matthew Dutson on 1/8/17.
//
// Contains a class used to generate random showers and perform a Monte Carlo simulation.

#ifndef monte_carlo_h
#define monte_carlo_h

#include <boost/property_tree/ptree.hpp>
#include <TRandom3.h>

#include "Simulator.h"

namespace cherenkov_library
{
    class MonteCarlo
    {
    public:

        /*
         * Takes a parsed XML object and attempts to extract required monte carlo parameters. If this fails, an
         * exception is thrown which specifies the name of the missing parameter.
         */
        MonteCarlo(boost::property_tree::ptree config_file);

        /*
         * Generates a random shower with a random direction, energy, and intensity profile. Allowed values and
         * distribution parameters are set in the configuration file.
         */
        Shower GenerateShower();

        /*
         * Constructs a shower object given a user-defined direction, impact parameter, and impact angle (the angle of
         * the point of closest approach).
         */
        Shower GenerateShower(TVector3 axis, double impact_param, double impact_angle);

        /*
         * Constructs a shower object given a user-defined direction, impact parameter, impact angle (the angle of the
         * point of closest approach), energy, and depth of first interaction.
         */
        Shower GenerateShower(TVector3 axis, double impact_param, double impact_angle, double energy);

    private:

        // The axis of the detector
        TVector3 detector_axis;

        // Parameters used to generate random showers in the Monte Carlo simulation
        TF1 energy_distribution;
        TF1 cosine_distribution;
        TF1 impact_distribution;
        double first_interact;
        double start_tracking;
        double n_max_ratio;

        // A general-purpose random number generator
        TRandom3 rng;

        // Atmospheric parameters
        double scale_height;
        double rho_0;
        double delta_0;

        // Parameters used when determining the depth of the shower maximum
        double x_max_1;
        double x_max_2;
        double x_max_3;
    };
}

#endif
