//
// Created by Matthew Dutson on 1/8/17.
//

#ifndef monte_carlo_h
#define monte_carlo_h

#include "simulator.h"
#include "boost/property_tree/ptree.hpp"
#include "TRandom3.h"

namespace cherenkov_library
{
    class MonteCarlo
    {
    public:

        void ParseFile(boost::property_tree::ptree config_file);

        Shower GenerateRandomShower();

    private:

        // Parameters used to generate random showers in the Monte Carlo simulation
        TF1 energy_distribution;
        TF1 cosine_distribution;
        TF1 impact_distribution;
        TF1 interact_distribution;
        double n_max_ratio;

        // A general-purpose random number generator
        TRandom3 rng;

        // Atmospheric parameters
        double scale_height;
        double rho_0;
        double delta_0;
    };
}

#endif
