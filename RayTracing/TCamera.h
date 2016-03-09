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
#include "TRawData.h"
#include "TSegmentedData.h"

class TCamera {
    
private:
    
    Double_t fHeight;
    
    Int_t fNumberTubesY;
    
    Double_t fWidth;
    
    Int_t fNumberTubesX;
    
    Double_t fPMTResponseTime;
    
    bool fTransparent;
    
    Double_t GetX(Int_t bin);
    
    Double_t GetY(Int_t bin);
    
    Int_t GetBin(Int_t x, Int_t y);
    
public:
    
    TCamera();
    
    TCamera(Double_t fHeight, Int_t numberTubesY, Double_t fWidth, Int_t numberTubesX, Double_t PMTResponseTime, bool transparent);
    
    TSegmentedData ParseData(TRawData data);
    
    void WriteDataToFile(TString filename, TSegmentedData);
    
    bool CheckCollision(TVector3 position);
    
};

#endif /* TCamera_h */
