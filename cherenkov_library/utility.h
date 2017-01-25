// common.h
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
// Contains a configuration manager and globally useful methods

#ifndef utility_h
#define utility_h

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "boost/property_tree/ptree.hpp"

namespace cherenkov_library
{

    /*
     * Reads the string and converts it to a TVector3 object.
     */
    TVector3 ToVector(std::string s);

    /*
     * A helper method which parses to double everything from the beginning of the string to the first occurence of the
     * specified character. The method then erases everything up to and including the specified character.
     */

    /*
     * Generates a randomly rotated vector perpendicular to the input. If the input vector is zero, (1, 0, 0) is
     * returned.
     */
    TVector3 RandomPerpendicularVector(TVector3 vec, TRandom3 rng);

    /*
     * Reads the file with the specified filename and parses it to XML. Throws exceptions with an informative message if
     * there is a problem reading or parsing the XML file
     */
    boost::property_tree::ptree ParseXMLFile(std::string filename);

    /*
     * Returns the speed of light measured in cm/s. Use this instead of the ROOT C() method.
     */
    double CentC();

    bool WithinXYDisk(TVector3 vec, double radius);

    /*
     * Writes the specified vector of arrays to a .csv file. Each array represents the contents of an individual row.
     */
    void WriteCSV(std::vector<std::vector<double>> data, std::vector<std::string> header, std::string filename);

    /*
     * Determines whether the xy projection of the vector lies within a disk centered at the origin.
     */
    bool WithinXYDisk(TVector3 vec, double radius);

    TRotation MakeRotation(double elevation_angle);
}


#endif
