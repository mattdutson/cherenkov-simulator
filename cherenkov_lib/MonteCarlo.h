// MonteCarlo.h
//
// Author: Matthew Dutson
//
// Definition of MonteCarlo class

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
    /*
     * A container for various user-defined XML parameters. Contains methods for performing the Monte Carlo simulation
     * and for generating random showers. Has a Simulator and Reconstructor as members, which are used during the Monte
     * Carlo.
     */
    class MonteCarlo
    {
    public:

        /*
         * Constructs the MonteCarlo by copying user-specified parameters from the parsed XML file.
         * TODO: Method should throw exceptions if parameters are out of range
         */
        explicit MonteCarlo(const boost::property_tree::ptree& config_file);

        /*
         * Performs the overall Monte Carlo simulation and writes results to a CSV file. A ROOT file is also written
         * which, for each shower, contains plots of the initial shower track, the post noise shower track, and the post
         * noise removal shower track.
         */
        void PerformMonteCarlo(std::string output_file) const;

        /*
         * Simulates and attempts reconstruction on a single shower, passed as a parameter. Writes various plots to the
         * current open file handle, and returns a Reconstructor::Result with reconstructed parameters. It is assumed
         * that a ROOT file will have been opened before calling this method.
         */
        Reconstructor::Result RunSingleShower(Shower shower, std::string ident) const;

        /*
         * Generates a Shower with a random position, direction, and energy. Allowed ranges of these parameters are
         * defined in the configuration file.
         */
        Shower GenerateShower() const;

        /*
         * Constructs a Shower given an axis direction, impact parameter, impact angle (angle of the point of closest
         * approach), and energy.
         */
        Shower GenerateShower(TVector3 axis, double im_par, double im_ang, double energy) const;

        /*
         * Parses the output file and configuration file from command line arguments, instantiates the MonteCarlo
         * object, and runs the PerformMonteCarlo method.
         */
        static int Run(int argc, const char* argv[]);

    private:

        friend class SampleEvents;

        int n_showers;
        double elevation;

        double energy_pow;
        double energy_min;
        double energy_max;
        double impact_min;
        double impact_max;
        double begn_depth;

        Simulator simulator;
        Reconstructor reconstructor;
    };
}

#endif
