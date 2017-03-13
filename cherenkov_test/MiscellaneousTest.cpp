// MiscellaneousTest.cpp
//
// Author: Matthew Dutson
//
// Contains miscellaneous tests.

#include <gtest/gtest.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Simulator.h"
#include "Utility.h"

using namespace std;
using namespace boost::property_tree;

namespace cherenkov_simulator
{
    /*
     * Reads a simple string to verify the behavior of the xml parser.
     */
    TEST(MiscellaneousTest, read_string)
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
    TEST(MiscellaneousTest, read_file)
    {
        ptree prop_tree;
        read_xml("SampleConfig.xml", prop_tree);
        ptree config = prop_tree.get_child("config");
        ASSERT_EQ(0.12, config.get<double>("surroundings.elevation_angle"));
    }
}