// Utility.h
//
// Author: Matthew Dutson
//
// Contains generally useful methods. Has no knowledge of Reconstructor, MonteCarlo, or Simulator.

#ifndef utility_h
#define utility_h

#include <boost/property_tree/ptree.hpp>
#include <string>
#include <TVector3.h>
#include <TRandom3.h>

namespace cherenkov_lib
{
    class Utility
    {
    public:

        /*
         * Reads the string and converts it to a TVector3 object.
         */
        static TVector3 ToVector(std::string s);

        /*
         * A helper method which parses to double everything from the beginning of the string to the first occurence of the
         * specified character. The method then erases everything up to and including the specified character.
         */
        static double ParseTo(std::string* s, char c);

        /*
         * Generates a randomly rotated vector perpendicular to the input. If the input vector is zero, (1, 0, 0) is
         * returned.
         */
        static TVector3 RandomPerpendicularVector(TVector3 vec, TRandom3* rng);

        /*
         * Reads the file with the specified filename and parses it to XML. Throws exceptions with an informative message if
         * there is a problem reading or parsing the XML file
         */
        static boost::property_tree::ptree ParseXMLFile(std::string filename);

        /*
         * Returns the speed of light measured in cm/s. Use this instead of the ROOT C() method.
         */
        static double CentC();

        /*
         * Writes the specified vector of arrays to a .csv file. Each array represents the contents of an individual row.
         */
        static void
        WriteCSV(std::vector<std::vector<double>> data, std::vector<std::string> header, std::string filename);

        /*
         * Determines whether the xy projection of the vector lies within a disk centered at the origin.
         */
        static bool WithinXYDisk(TVector3 vec, double radius);

        static TRotation MakeRotation(double elevation_angle);

        static bool Above(TVector3 reference, TVector3 other);

        static double RandLinear(TRandom3* rng, double max);
    };
}


#endif
