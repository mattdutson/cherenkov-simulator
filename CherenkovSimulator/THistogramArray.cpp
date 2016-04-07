//
//  TPixelHistograms.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/29/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "THistogramArray.h"

THistogramArray::THistogramArray() {
    fHistograms = std::list<TPixelData>();
    fNBins = 0;
}

void THistogramArray::AddHistogram(Double_t x, Double_t y, TH1D histogram) {
    TPixelData pixelData = TPixelData(x, y, histogram);
    fHistograms.push_back(pixelData);
}

Int_t THistogramArray::GetNBins() {
    return fNBins;
}

std::list<TPixelData>::iterator THistogramArray::Begin() {
    return fHistograms.begin();
}

std::list<TPixelData>::iterator THistogramArray::End() {
    return fHistograms.end();
}

void THistogramArray::WriteToFile(TString filename) {
    TFile file(filename, "RECREATE");
    std::list<TPixelData>::iterator iter;
    for (iter = fHistograms.begin(); iter != fHistograms.end(); iter++) {
        (*iter).Write();
    }
    file.Close();
}