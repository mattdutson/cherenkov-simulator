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
#include "TVector3.h"

class TCamera {
    
private:
    
    Double_t fHeight;
    
    Int_t fNumberTubesY;
    
    Double_t fWidth;
    
    Int_t fNumberTubesX;
    
    Double_t fPMTResponseTime;
    
    Double_t fMinTime;
    
    Double_t fMaxTime;
    
    Double_t GetPMTX(Int_t xIndex);
    
    Double_t GetPMTY(Int_t yIndex);
    
    Double_t GetXBin(Double_t x);
    
    Double_t GetYBin(Double_t y);
    
public:
    
    TCamera();
    
    TCamera(Double_t fHeight, Int_t numberTubesY, Double_t fWidth, Int_t numberTubesX, Double_t PMTResponseTime);
    
    std::vector<Double_t>*** ParseData(std::vector<Double_t> x, std::vector<Double_t> y, std::vector<Double_t> time);
    
    void WriteDataToFile(TString filename, std::vector<Double_t>*** data);
    
    bool CheckCollision(TVector3 position);
    
};

#endif /* TCamera_h */
