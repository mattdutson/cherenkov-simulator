//
//  TCamera.hpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/23/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TCamera_h
#define TCamera_h

#include <stdio.h>
#include "TMath.h"

class TCamera {
    
private:
    
    Double_t fHeight;
    
    Int_t fNumberTubesY;
    
    Double_t fWidth;
    
    Int_t fNumberTubesX;
    
    Double_t fPMTResponseTime;
    
    Double_t fMinTime;
    
    Double_t fMaxTime;
    
public:
    
    TCamera();
    
    TCamera(Double_t fHeight, Int_t numberTubesY, Double_t fWidth, Int_t numberTubesX, Double_t PMTResponseTime);
    
    std::vector<Double_t>*** ParseData(std::vector<Double_t> x, std::vector<Double_t> y, std::vector<Double_t> time);
    
    void WriteDataToFile(TString filename, std::vector<Double_t>*** data);
    
};

#endif /* TCamera_h */
