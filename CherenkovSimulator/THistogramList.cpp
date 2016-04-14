/*
 * CherenkovSimulator - TCoordinates.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of THistogramList.
 */

#include "THistogramList.hpp"

THistogramList::THistogramList() {
    fHistograms = std::list<TPixelData>();
    fNBins = 0;
}

void THistogramList::AddHistogram(Double_t x, Double_t y, TH1D histogram) {
    TPixelData pixelData = TPixelData(x, y, histogram);
    fHistograms.push_back(pixelData);
}

Int_t THistogramList::GetNBins() {
    return fNBins;
}

std::list<TPixelData>::iterator THistogramList::Begin() {
    return fHistograms.begin();
}

std::list<TPixelData>::iterator THistogramList::End() {
    return fHistograms.end();
}

void THistogramList::WriteToFile(TString filename) {
    TFile file(filename, "RECREATE");
    std::list<TPixelData>::iterator iter;
    for (iter = fHistograms.begin(); iter != fHistograms.end(); iter++) {
        (*iter).Write();
    }
    file.Close();
}