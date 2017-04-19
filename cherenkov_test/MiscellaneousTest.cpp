// MiscellaneousTest.cpp
//
// Author: Matthew Dutson
//
// Contains miscellaneous tests.

#include <gtest/gtest.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <TFile.h>
#include <TH1I.h>

#include "MonteCarlo.h"

using namespace std;
using namespace boost::property_tree;

namespace cherenkov_simulator
{
    /*
     * Reads a simple string to verify the behavior of the xml parser.
     */
    TEST(MiscellaneousTest, ReadString)
    {
        stringstream file;
        ptree prop_tree;
        file << "<configuration><property>10</property></configuration>";
        read_xml(file, prop_tree);
        ptree config = prop_tree.get_child("configuration");
        ASSERT_EQ(10, config.get<int>("property"));
    }

    /*
     * Reads a simple XML file from a filename to verify the behavior of the ptree xml parser.
     */
    TEST(MiscellaneousTest, ReadFile)
    {
        ptree prop_tree;
        read_xml("SampleConfig.xml", prop_tree);
        ptree config = prop_tree.get_child("config");
        ASSERT_EQ(0.12, config.get<double>("surroundings.elevation_angle"));
    }

    TEST(MiscellaneousTest, MonteCarlo)
    {
        const char* argv[] = {"", "Output", "../../../Config.xml", "1"};
        cherenkov_simulator::MonteCarlo::Run(4, argv);
    }

    TEST(MiscellaneousTest, PowerLaw)
    {
        TFile file("PowerLaw.root", "RECREATE");
        TH1I power_histo = TH1I("power_law", "Power Law Value", 10000, 1e2, 1e8);
        for (int i = 0; i < 1e8; i++)
        {
            double value = Utility::RandPower(1e2, 1e8, -1);
            power_histo.Fill(value);
        }
        power_histo.Write("power_histo");
    }
}