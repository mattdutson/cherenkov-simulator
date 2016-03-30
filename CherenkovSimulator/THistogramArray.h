//
//  TPixelHistograms.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/29/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef THistogramArray_h
#define THistogramArray_h

#include "TH1D.h"
#include "TSegmentedData.h"
#include <stdio.h>

class THistogramArray {
private:
    
    Int_t fNBins;
    
    Double_t fMinTime;
    
    Double_t fMaxTime;
    
    TH1D** fHistograms;
    
public:
    
    THistogramArray(Int_t nBins);
    
    void SetHistogram(Int_t bin, TH1D histogram);
    
    TH1D GetHistogram(Int_t bin);
    
    Int_t GetNBins();
    
};

#endif /* THistogramArray.h */
