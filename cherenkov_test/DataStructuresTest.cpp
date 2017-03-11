// DataStructuresTest.cpp
//
// Author: Matthew Dutson
//
// Tests of DataStructures.h

#include <gtest/gtest.h>
#include <TFile.h>

#include "DataStructures.h"
#include "Analysis.h"

using namespace cherenkov_lib;

namespace cherenkov_test
{
    /*
     * Outputs a 2D histogram to a file. This histogram will contain a map of valid pixels.
     */
    TEST(DataStructuresTest, MakeValidMap)
    {
        PhotonCount counter = PhotonCount(200, 0.0, 10.0e-9, 0.007, 1.0);
        TH2C map = Analysis::GetValidMap(counter);
        TFile file("pixel_map.root", "RECREATE");
        map.Write("valid");
    }
}
