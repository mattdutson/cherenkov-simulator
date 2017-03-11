//
// Created by Matthew Dutson on 1/26/17.
//

#include <gtest/gtest.h>
#include <TFile.h>

#include "DataStructures.h"
#include "Helper.h"
#include "Analysis.h"

using namespace cherenkov_library;

namespace cherenkov_tests
{
    /*
     * Outputs a 2D histogram to a file. This histogram will contain a map of valid pixels.
     */
    TEST(DataContainersTest, MakeValidMap)
    {
    PhotonCount counter = PhotonCount(200, 0.0, 10.0e-9, 0.007, 1.0);
        TH2C map = DataAnalysis::GetValidMap(counter);
    TFile file("../../cherenkov_test/pixel_map.root", "RECREATE");
        map.Write("valid");
    }
}
