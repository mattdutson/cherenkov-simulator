//
// Created by Matthew Dutson on 10/18/16.
//

#include <gtest/gtest.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "simulator.h"

using namespace std;
using namespace boost::property_tree;

/*
 * Reads a simple XML string to verify the behavior of the ptree xml parser.
 */
TEST(config_file, read_string)
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
TEST(config_file, read_file)
{
    ptree prop_tree;
    read_xml("../../cherenkov_tests/sample_config.xml", prop_tree);
    ptree config = prop_tree.get_child("configuration");
    ASSERT_EQ(0.125, config.get<double>("elevation_angle"));
}