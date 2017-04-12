// MonteCarlo.h
//
// Author: Matthew Dutson
//
// Contains a class used to generate random showers for the Monte Carlo simulation.

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <boost/property_tree/ptree.hpp>
#include <TF1.h>
#include <TRandom3.h>

#include "Geometric.h"
#include "Simulator.h"
#include "Reconstructor.h"
#include "Utility.h"

namespace cherenkov_simulator
{
    class MonteCarlo
    {
    public:

        /*
         * Constructs the MonteCarlo from values in the configuration tree.
         */
        MonteCarlo(const boost::property_tree::ptree& config_file);

        /*
         * Performs the overall Monte Carlo simulation and writes results to the console.
         */
        void PerformMonteCarlo(std::string out_file);

        /*
         * Generates a random shower with a random direction, energy, and intensity profile. Allowed values and
         * distribution parameters are set in the configuration file.
         */
        Shower GenerateShower();

        /*
         * Constructs a shower object given a user-defined direction, impact parameter, impact angle (the angle of the
         * point of closest approach), energy, and depth of first interaction.
         */
        Shower GenerateShower(TVector3 axis, double impact_param, double impact_angle, double energy);

        /*
         * Parses the output file and configuration file from command line arguments, instantiates the MonteCarlo
         * object, and runs PerformMonteCarlo method.
         */
        static int Run(int argc, const char* argv[]);

    private:

        // Used to calculate N_max from the energy - eV
        constexpr static double n_max_ratio = 1.39e9;

        // Parameters used when determining the depth of the shower maximum - cgs
        constexpr static double x_max_1 = 725.0;
        constexpr static double x_max_2 = 55.0;
        constexpr static double x_max_3 = 18.0;

        // Atmospheric parameters - cgs
        constexpr static double scale_height = 841300;
        constexpr static double rho_sea = 0.001225;
        constexpr static double refrac_sea = 1.00029;

        // Non-constant atmospheric parameters (depend on the elevation) - cgs
        double rho_0;
        double delta_0;

        // Parameters used to generate random showers in the Monte Carlo simulation - eV, cgs
        TF1 energy_distribution;
        TF1 cosine_distribution;
        TF1 impact_distribution;
        double start_tracking;

        // The number of showers to simulate when PerformMonteCarlo is called
        int n_showers;

        // A general-purpose random number generator
        TRandom3 rng;

        // The Simulator and Reconstructor objects used when PerformMonteCarlo is called
        Simulator simulator;
        Reconstructor reconstructor;
    };
}

#endif
