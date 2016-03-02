//
//  TDataCollection.hpp
//  RayTracing
//
//  Created by Matthew Dutson on 3/1/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TDataCollection_h
#define TDataCollection_h

#include <stdio.h>
#include "TMath.h"

class TDataCollection {
    
private:
    
    std::vector<Double_t> fXData;
    
    std::vector<Double_t> fYData;
    
    std::vector<Double_t> fTData;
    
    Double_t fXMin;
    
    Double_t fXMax;
    
    Double_t fYMin;
    
    Double_t fYMax;
    
    Double_t fTMin;
    
    Double_t fTMax;
    
public:
    
    TDataCollection();
    
    void PushBack(Double_t x, Double_t y, Double_t t);
    
    void Clear();
    
    std::vector<Double_t> GetXData();
    
    std::vector<Double_t> GetYData();
    
    std::vector<Double_t> GetTData();
    
};

#endif /* TDataCollection_h */
