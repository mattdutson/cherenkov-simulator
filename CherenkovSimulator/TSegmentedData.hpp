/*
 * CherenkovSimulator - TSegmentedData.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 *
 */

#ifndef TSegmentedData_h
#define TSegmentedData_h

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
