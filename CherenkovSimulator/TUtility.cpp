/*
 * CherenkovSimulator - TRay.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 *
 */

#include "TUtility.h"

Double_t TUtility::SumArray(std::vector<Double_t> array) {
    Double_t sum = 0;
    for (Double_t d: array) {
        sum += d;
    }
    return sum;
}

void TUtility::WriteHistogramFile(TString filename, THistogramArray histograms) {
    TFile file(filename, "RECREATE");
    for (std::list<TPixelData>::iterator iter = histograms.Begin(); iter != histograms.End(); iter++) {
        (*iter).Write();
    }
    file.Close();
}