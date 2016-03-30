//
//  TPixelHistograms.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/29/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TPixelHistograms.h"

TPixelHistograms::TPixelHistograms(Int_t nBins) {
    fHistograms = new TH1D*[nBins];
    for (Int_t i = 0; i < nBins; i++) {
        fHistograms[i] = new TH1D();
    }
    fNBins = nBins;
    fMinTime = 1e100;
    fMaxTime = -1e100;
}

void TPixelHistograms::SetHistogram(Int_t bin, TH1D histogram) {
    fHistograms[bin] = &histogram;
}