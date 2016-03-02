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
    
public:
    
    TDataCollection();
    
    void PushBack(Double_t x, Double_t y, Double_t t);
    
    void Clear();
    
    std::vector<Double_t> GetXData();
    
    std::vector<Double_t> GetYData();
    
    std::vector<Double_t> GetTData();
    
    Double_t GetX(Int_t index);
    
    Double_t GetY(Int_t index);
    
    Double_t GetT(Int_t index);
    
    Double_t Size();
    
};

#endif /* TDataCollection_h */
