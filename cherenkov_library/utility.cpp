// common.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "utility.h"
#include "TMath.h"

using namespace boost::program_options;
using namespace std;
using boost::property_tree::ptree;

namespace cherenkov_simulator
{
    TVector3 ToVector(string s)
    {
        TVector3 output = TVector3();
        int comma = s.find(',');
        output.SetX(stod(s.substr(1, comma)));
        s.erase(0, comma);
        comma = s.find(',');
        output.SetY(stod(s.substr(0, comma)));
        s.erase(0, comma);
        output.SetZ(stod(s.substr(0, s.size() - 1)));
        return output;
    }

    TVector3 RandomPerpendicularVector(TVector3 vec, TRandom3 rng)
    {
        if (vec.X() == 0 && vec.Y() == 0 && vec.Z() == 0)
        {
            return TVector3(1, 0, 0);
        } else
        {
            TVector3 other_vec = vec + TVector3(1, 0, 0);
            TVector3 normal = (vec.Cross(other_vec)).Unit();
            normal.Rotate(rng.Uniform(2 * TMath::Pi()), vec);
            return normal;
        }
    }

    ptree ParseXMLFile(string filename)
    {
        // Try opening the specified file
        ifstream config_file = ifstream();
        try
        {
            config_file.open(filename);
        }
        catch (...)
        {
            cout << "The file " << filename << " could not be opened. Check the path." << endl;
            return -1;
        }

        // Try to construct the simulator using parameters in the configuration file
        Simulator sim = Simulator();
        try
        {
            sim.ParseFile(config_file);
        }
        catch (runtime_error e)
        {
            cout << e.what() << endl;
            return -1;
        }
        catch (...)
        {
            cout << "An unknown exception occured." << endl;
            return -1;
        }

        // Parse the file to XML.
        ptree xml_file = ptree();
        try
        {
            read_xml(config_file, xml_file);
        }
        catch (...)
        {
            throw std::runtime_error("There was a problem parsing the file to XML. Check for syntax errors.");
        }

        bool Utility::WithinXYDisk(TVector3 vec, double radius)
        {
            double xy_radius = Sqrt(vec.X() * vec.X() + vec.Y() * vec.Y());
            if (xy_radius < radius)
            {
                return true;
            } else
            {
                return false;
            }
        }
    }

    double CentC()
    {
        return TMath::C() * 100.0;
    }
}
