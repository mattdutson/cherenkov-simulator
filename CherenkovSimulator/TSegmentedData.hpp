/*
 * CherenkovSimulator - TSegmentedData.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains raw data which has been separated into bins according to the locations of photomultipler tubes.
 */

#ifndef TSegmentedData_hpp
#define TSegmentedData_hpp

#include <stdio.h>
#include "TVector3.h"

class TSegmentedData {
    
private:
    
    Int_t fNBins;
    
    Double_t fMinTime;
    
    Double_t fMaxTime;
    
    std::vector<Double_t>** fSegmentedData;

public:
    
    TSegmentedData(Int_t nBins);
    
    void AddPoint(Double_t time, Int_t bin);
    
    Double_t GetMaxTime();
    
    Double_t GetMinTime();
    
    std::vector<Double_t>* GetSegment(Int_t bin);
    
    Int_t GetNBins();
};

#endif /* TSegmentedData_h */
