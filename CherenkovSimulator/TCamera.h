//
//  TCamera.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TCamera_h
#define TCamera_h

#include "TPlane3.h"
#include "TSegmentedData.h"
#include "TRawData.h"
#include "TRay.h"

class TCamera {
    
private:
    
    Double_t fFocalLength;
    
    Double_t fHeight;
    
    Int_t fNumberTubesY;
    
    Double_t fWidth;
    
    Int_t fNumberTubesX;
    
    Double_t fPMTResolution;
    
    Bool_t fCheckBackCollision = false;
    
    Int_t GetBin(TVector3 position);
    
public:
    
    TCamera();
    
    TCamera(Double_t focalLength, Double_t width, Int_t numberTubesX, Double_t height, Int_t numberTubesY, Double_t PMTResolution, Bool_t checkBackCollision);
    
    Double_t FocalLength();
    
    Bool_t CheckCollision(TVector3 position);
    
    TSegmentedData ParseData(TRawData data);
    
    TVector3 GetViewDirection(Int_t bin);
};

#endif /* TCamera_h */
