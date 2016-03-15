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
#include "TAnalysis.h"

class TCamera {
    
private:
    
    Double_t fFocalLength;
    
    TPlane3 fFocalPlane;
    
    Double_t fHeight;
    
    Int_t fNumberTubesY;
    
    Double_t fWidth;
    
    Int_t fNumberTubesX;
    
    Double_t fPMTResponseTime;
    
    Int_t GetBin(TVector3 position);
    
public:
    
    TPlane3 FocalPlane();
    
    Bool_t CheckCollision(TVector3 position);
    
    TSegmentedData ParseData(TRawData data);
    
    TVector3 GetViewDirection(Int_t bin);
};

#endif /* TCamera_hpp */
