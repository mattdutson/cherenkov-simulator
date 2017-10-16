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

#include "Utility.h"
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
        void PerformMonteCarlo(std::string out_file) const;

        /*
         * Generates a random shower with a random direction, energy, and intensity profile. Allowed values and
         * distribution parameters are set in the configuration file.
         */
        Shower GenerateShower() const;

        /*
         * Constructs a shower object given a user-defined direction, impact parameter, impact angle (the angle of the
         * point of closest approach), energy, and depth of first interaction.
         */
        Shower GenerateShower(TVector3 axis, double impact_param, double impact_angle, double energy) const;

        /*
         * Parses the output file and configuration file from command line arguments, instantiates the MonteCarlo
         * object, and runs PerformMonteCarlo method.
         */
        static int Run(int argc, const char* argv[]);

    private:

        // Parameters used to generate random showers in the Monte Carlo simulation - eV, cgs
        double elevation;
        double energy_pow;
        double e_min;
        double e_max;
        double impact_min;
        double impact_max;
        double start_tracking;

        // The number of showers to simulate when PerformMonteCarlo is called
        int n_showers;

        // The Simulator and Reconstructor objects used when PerformMonteCarlo is called
        Simulator simulator;
        Reconstructor reconstructor;
    };
}

#endif
