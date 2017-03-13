// DataStructuresTest.cpp
//
// Author: Matthew Dutson
//
// Tests of DataStructures.h

#include <gtest/gtest.h>
#include <TFile.h>

#include "DataStructures.h"
#include "Analysis.h"

namespace cherenkov_simulator
{
    /*
     * Note: this class will be able to access private member of the PhotonCount and Iterator classes
     */
    class DataStructuresTest : public ::testing::Test
    {
    };

    /*
     * Outputs a 2D histogram to a file. This histogram will contain a map of valid pixels.
     */
    TEST_F(DataStructuresTest, MakeValidMap)
    {
        PhotonCount::Params params;
        params.n_pixels = 200;
        params.bin_size = 10.0e-9;
        params.angular_size = 0.007;
        params.linear_size = 1.0;
        PhotonCount counter = PhotonCount(params, 0.0);
        TH2C map = Analysis::GetValidMap(counter);
        TFile file("pixel_map.root", "RECREATE");
        map.Write("valid");
    }

}
