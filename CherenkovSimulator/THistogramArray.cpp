//
//  TPixelHistograms.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/29/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "THistogramArray.h"

THistogramArray::THistogramArray(Int_t nBins) {
    fHistograms = new TH1*[nBins];
    fNBins = nBins;
    fMinTime = 1e100;
    fMaxTime = -1e100;
}

void THistogramArray::SetHistogram(Int_t bin, TH1* histogram) {
    fHistograms[bin] = histogram;
}

Int_t THistogramArray::GetNBins() {
    return fNBins;
}

TH1* THistogramArray::GetHistogram(Int_t bin) {
    return fHistograms[bin];
}