//
//  TPixelHistograms.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/29/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TPixelHistograms_h
#define TPixelHistograms_h

#include "TH1D.h"
#include "TSegmentedData.h"
#include <stdio.h>

class TPixelHistograms {
private:
    
    Int_t fNBins;
    
    Double_t fMinTime;
    
    Double_t fMaxTime;
    
    TH1D** fHistograms;
    
public:
    
    TPixelHistograms(Int_t nBins);
    
    void SetHistogram(Int_t bin, TH1D histogram);
    
};

#endif /* TPixelHistograms_h */
