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
#include "THistogramArray.h"
#include "TRawData.h"
#include "TRay.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

class TCamera {
    
private:
    
    Double_t fFocalLength;
    
    Double_t fHeight;
    
    Int_t fNumberTubesY;
    
    Double_t fWidth;
    
    Int_t fNumberTubesX;
    
    Double_t fPMTResolution;
    
    TF1 fResponseFunction;
    
    Bool_t fCheckBackCollision = false;
    
    Int_t GetBin(TVector3 position);
    
public:
    
    TCamera();
    
    TCamera(Double_t focalLength, Double_t width, Int_t numberTubesX, Double_t height, Int_t numberTubesY, Double_t PMTResolution, TF1 responseFunction, Bool_t checkBackCollision);
    
    Double_t FocalLength();
    
    Bool_t CheckCollision(TVector3 position);
    
    TSegmentedData SegmentedData(TRawData data);
    
    THistogramArray PixelHistograms(TSegmentedData data);
    
    THistogramArray VoltageOutput(THistogramArray histograms);
    
    TVector3 GetViewDirection(Int_t bin);
    
    TVector2 GetPixelPosition(Int_t bin);
    
    void WriteDataToFile(TString filename, TSegmentedData parsedData);
};

#endif /* TCamera_h */
