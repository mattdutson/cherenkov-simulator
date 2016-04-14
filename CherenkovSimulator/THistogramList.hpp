/*
 * CherenkovSimulator - THistogramList.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains a linked list of histograms which is intented for traversal in order. This is mostly used for output voltage data.
 */

#ifndef THistogramList_hpp
#define THistogramList_hpp

#include "TSegmentedData.hpp"
#include "TPixelData.hpp"

#include "TH1D.h"
#include "TFile.h"
#include <stdio.h>
#include <list>

class THistogramList {
    
private:
    
    Int_t fNBins;
    
    std::list<TPixelData> fHistograms;

public:
    
    THistogramList();
    
    void AddHistogram(Double_t x, Double_t y, TH1D histogram);
    
    std::list<TPixelData>::iterator Begin();
    
    std::list<TPixelData>::iterator End();
    
    Int_t GetNBins();
    
    void WriteToFile(TString filename);
};

#endif
