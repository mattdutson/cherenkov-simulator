// Utility.h
//
// Author: Matthew Dutson
//
// Contains generally useful methods. Has no knowledge of Reconstructor, MonteCarlo, or Simulator.

#ifndef UTILITY_H
#define UTILITY_H

#include <boost/property_tree/ptree.hpp>
#include <string>
#include <vector>
#include <TVector3.h>

namespace cherenkov_simulator
{
    // Physics constants - MeV, cgs
    const double mass_e = 0.511;
    const double fine_struct = 0.007297;
    const double c_cent = 2.998e10;

    // Fixed parameters in the Gaisser-Hillas profile - cgs.
    const double x_0 = -70.0;
    const double gh_lambda = 70.0;

    // Used to calculate N_max from the energy - eV
    const double n_max_ratio = 1.39e9;

    // Parameters used when determining the depth of the shower maximum - cgs
    const double x_max_1 = 725.0;
    const double x_max_2 = 55.0;
    const double x_max_3 = 18.0;

    // Atmospheric parameters - cgs
    const double atm_h = 841300;
    const double rho_sea = 0.001225;
    const double refrac_sea = 1.00029;

    // Levels of night sky background noise - cgs, sr
    const double global_sky_noise = 4.924e6;
    const double global_ground_noise = 4.924e5;

    // Parameters in the Cherenkov yield - cgs
    const double lambda_min = 3.0e-5;
    const double lambda_max = 4.0e-5;

    // Parameters in the electron energy spectrum - MeV
    const double fe_a11 = 6.42522;
    const double fe_a12 = 1.53183;
    const double fe_a21 = 168.168;
    const double fe_a22 = 42.1368;
    const double fe_k0 = 1.45098e-1;
    const double fe_k1 = 6.20114;
    const double fe_k2 = -5.96851e-1;

    // Miscellaneous parameters - K
    const double atm_temp = 273.0;
    const double refrac_lens = 1.52;

    // Parameters in the fluorescence yield - MeV, cgs, K
    const double fluor_a1 = 890.0;
    const double fluor_a2 = 550.0;
    const double fluor_b1 = 1850.0;
    const double fluor_b2 = 6500.0;
    const double dep_1_4 = 1.6;

    // Parameters in the effective ionization loss rate - MeV, cgs
    const double ion_c1 = 3.90883;
    const double ion_c2 = 1.05301;
    const double ion_c3 = 9.91717;
    const double ion_c4 = 2.41715;
    const double ion_c5 = 0.13180;

    // Parameters used when calculating theta_c in the Cherenkov angular distribution - MeV
    const double ckv_k1 = 0.83;
    const double ckv_k2 = -0.67;

    // Parameters which describe inefficiencies in the equipment
    const double mirror_reflect = 0.80; // 1.0; // 0.80
    const double filter_transmit = 1.0;
    const double quantum_eff = 0.15; // 1.0; // 0.15
    
    
    typedef std::vector<bool> Bool1D;
    typedef std::vector<std::vector<bool>> Bool2D;
    typedef std::vector<std::vector<std::vector<bool>>> Bool3D;

    typedef std::vector<short> Short1D;
    typedef std::vector<std::vector<short>> Short2D;
    typedef std::vector<std::vector<std::vector<short>>> Short3D;

    typedef std::vector<double> Double1D;

    class Utility
    {
    public:
        
        /*
         * Reads the string and converts it to a TVector3 object.
         */
        static TVector3 ToVector(std::string s);

        /*
         * Reads the file with the specified filename and parses it to XML. Throws exceptions with an informative
         * message if there is a problem reading or parsing the XML file
         */
        static boost::property_tree::ptree ParseXMLFile(std::string filename);

        /*
         * Determines whether the xy projection of the vector lies within a disk centered at the origin.
         */
        static bool WithinXYDisk(TVector3 vec, double radius);

        /*
         * Constructs the rotation used by MonteCarlo and Simulator classes.
         */
        static TRotation MakeRotation(double elevation_angle);

        /*
         * Generates a randomly rotated vector perpendicular to the input. If the input vector is zero, (1, 0, 0) is
         * returned.
         */
        static TVector3 RandNormal(TVector3 vec);

        /*
         * Returns a random, linearly distributed value constrained between zero and some maximum.
         */
        static double RandLinear(double min, double max);

        /*
         * Generates a random angle on (0, pi) weighted by a cosine.
         */
        static double RandCosine();

        /*
         * Generates a random number according to a power law distribution.
         */
        static double RandPower(double min, double max, double pow);

        /*
         * Returns an integer which is randomly rounded up or down from the input double based on its decimal. For
         * instance, 3.2 would be rounded up to 4 20% of the time and down to 3 80% of the time.
         */
        static int RandomRound(double value);

        /*
         * Calculates the percent error between the actual and expected values. If the expected value is zero, the
         * actual value is returned.
         */
        static double PercentError(double actual, double expected);

        /*
         * Converts a centimeter double into a kilometer string.
         */
        static std::string KmString(double cent);

    private:

        /*
         * A helper method which parses to double everything from the beginning of the string to the first occurence of the
         * secified character. The method then erases everything up to and including the specified character.
         */
        static double ParseTo(std::string& s, char c);
    };
}


#endif
