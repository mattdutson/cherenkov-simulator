//
// Created by Matthew Dutson on 1/8/17.
//

#ifndef monte_carlo_h
#define monte_carlo_h

#include "simulator.h"
#include "boost/property_tree/ptree.hpp"

namespace cherenkov_simulator
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
    };
}

#endif
