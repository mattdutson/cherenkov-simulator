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
#include "TDataCollection.h"

class TCamera {
    
private:
    
    Double_t fHeight;
    
    Int_t fNumberTubesY;
    
    Double_t fWidth;
    
    Int_t fNumberTubesX;
    
    Double_t fPMTResponseTime;
    
    bool fTransparent;
    
    Double_t fMinTime;
    
    Double_t fMaxTime;
    
    Double_t GetPMTX(Int_t xIndex);
    
    Double_t GetPMTY(Int_t yIndex);
    
    Double_t GetXBin(Double_t x);
    
    Double_t GetYBin(Double_t y);
    
public:
    
    TCamera();
    
    TCamera(Double_t fHeight, Int_t numberTubesY, Double_t fWidth, Int_t numberTubesX, Double_t PMTResponseTime, bool transparent);
    
    std::vector<Double_t>*** ParseData(TDataCollection data);
    
    void WriteDataToFile(TString filename, std::vector<Double_t>*** parsedData);
    
    bool CheckCollision(TVector3 position);
    
};

#endif /* TCamera_h */
