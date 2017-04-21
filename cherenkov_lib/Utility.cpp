// Utility.cpp
//
// Author: Matthew Dutson
//
// Implementation of Utility.h

#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <TRotation.h>
#include <TRandom3.h>

#include "Utility.h"

using namespace std;
using namespace TMath;

using boost::property_tree::ptree;

namespace cherenkov_simulator
{
    TVector3 Utility::ToVector(string s)
    {
        // Clear out everything before the first parenthesis.
        size_t current = s.find('(');
        s.erase(0, current + 1);

        // Pull double values from the vector.
        TVector3 output = TVector3();
        output.SetX(ParseTo(s, ','));
        output.SetY(ParseTo(s, ','));
        output.SetZ(ParseTo(s, ')'));
        return output;
    }

    ptree Utility::ParseXMLFile(string filename)
    {
        // Try opening the specified file
        ifstream config_file = ifstream(filename);
        if (config_file.fail())
        {
            std::string message = "The file " + filename + " could not be opened. Check the path.";
            throw std::runtime_error(message);
        }

        // Parse the file to XML.
        try
        {
            ptree xml_file = ptree();
            read_xml(config_file, xml_file);
            return xml_file;
        }
        catch (boost::property_tree::xml_parser_error)
        {
            throw std::runtime_error("There was a problem parsing the file to XML. Check for syntax errors.");
        }
    }

    TRotation Utility::MakeRotation(double elevation_angle)
    {
        TRotation rotate = TRotation();
        rotate.RotateX(-PiOver2() + elevation_angle);
        return rotate;
    }

    bool Utility::WithinXYDisk(TVector3 vec, double radius)
    {
        return Sqrt(Sq(vec.X()) + Sq(vec.Y())) < radius;
    }

    TVector3 Utility::RandNormal(TVector3 vec)
    {
        if (vec.Mag2() == 0)
        {
            return TVector3(1, 0, 0);
        }
        else
        {
            TVector3 other_vec = vec + TVector3(1, 0, 0);
            TVector3 normal = (vec.Cross(other_vec)).Unit();
            normal.Rotate(gRandom->Uniform(2 * TMath::Pi()), vec);
            return normal;
        }
    }

    double Utility::RandLinear(double min, double max)
    {
        if (min < 0 || max < 0) throw std::runtime_error("The bounds must be non-negative");
        if (min >= max) throw std::runtime_error("The min bound must be less than the max bound");
        return Sqrt((Sq(max) - Sq(min)) * gRandom->Rndm() + Sq(min));
    }

    double Utility::RandCosine()
    {
        return ASin(gRandom->Rndm());
    }

    double Utility::RandPower(double min, double max, double pow)
    {
        if (min <= 0 || max <= 0) throw std::runtime_error("The bounds must be positive");
        if (min >= max) throw std::runtime_error("The min bound must be less than the max bound");
        if (pow == -1)
        {
            return min * Power(max / min, gRandom->Rndm());
        }
        else
        {
            double a = Power(min, pow + 1);
            double b = Power(max, pow + 1);
            return Power((b - a) * gRandom->Rndm() + a, 1.0 / (pow + 1));
        }
    }

    int Utility::RandomRound(double value)
    {
        double decimal = value - Floor(value);
        int base = (int) (value - decimal);
        if (gRandom->Rndm() < decimal) return base + 1;
        else return base;
    }

    double Utility::ParseTo(string& s, char c)
    {
        size_t index = s.find(c);
        double out = stod(s.substr(0, index));
        s.erase(0, index + 1);
        return out;
    }

    double Utility::PercentError(double actual, double expected)
    {
        return expected == 0 ? Abs(actual) : Abs((actual - expected) / expected);
    }

    string Utility::KmString(double cent)
    {
        return to_string(cent / 1e5);
    }
}
