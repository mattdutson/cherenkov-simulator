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

namespace cherenkov_test
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
read_xml("../../cherenkov_test/sample_config.xml", prop_tree);
        ptree config = prop_tree.get_child("configuration");
        ASSERT_EQ(0.125, config.get<double>("elevation_angle"));
    }

    /*
     * Writes a simple CSV from a vector of arrays of strings.
     */
    TEST(MiscellaneousTest, write_string_csv)
    {
        vector<string> header = {"ID", "ISBN", "Price"};
        vector<double> row1 = {1, 1234567890, 19.99};
        vector<double> row2 = {2, 9876543210, 8.72};
        vector<double> row3 = {3, 2222222222, 32.95};
        vector<vector<double>> all_rows = vector<vector<double>>();
        all_rows.push_back(row1);
        all_rows.push_back(row2);
        all_rows.push_back(row3);
        cherenkov_lib::Utility::WriteCSV(all_rows, header,
                                         "../../cherenkov_test/sample_csv.csv");
    }
}