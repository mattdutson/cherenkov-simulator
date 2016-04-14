/*
 * CherenkovSimulator - THistogram.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 *
 */

#ifndef THistogramArray_h
#define THistogramArray_h

#include "TH1D.h"
#include "TSegmentedData.h"
#include "TFile.h"
#include "TPixelData.h"
#include <stdio.h>
#include <list>

class THistogramArray {
    
private:
    
    Int_t fNBins;
    
    std::list<TPixelData> fHistograms;

public:
    
    THistogramArray();
    
    void AddHistogram(Double_t x, Double_t y, TH1D histogram);
    
    std::list<TPixelData>::iterator Begin();
    
    std::list<TPixelData>::iterator End();
    
    Int_t GetNBins();
    
    void WriteToFile(TString filename);
};

#endif /* THistogramArray.h */
