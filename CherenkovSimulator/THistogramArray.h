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
#include "TFile.h"
#include <stdio.h>
#include <list>

class TPixelData: public TH1D {
    
private:
    
    Double_t fX;
    
    Double_t fY;
    
public:
    
    void SetPosition(Double_t x, Double_t y);
};

class THistogramArray {
    
private:
    
    Int_t fNBins;
    
    std::list<TPixelData> fHistograms;

public:
    
    THistogramArray();
    
    void AddHistogram(Double_t x, Double_t y, TH1D histogram);
    
    std::list<TPixelData>::iterator GetHistogram(Int_t bin);
    
    Int_t GetNBins();
    
    void WriteToFile(TString filename);
};

#endif /* THistogramArray.h */
