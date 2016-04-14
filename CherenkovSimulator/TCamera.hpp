/*
 * CherenkovSimulator - TCamera.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Represents the detection apparatus of the cosmic ray observatory. The detection apparatus consists of a rectangular array of photomultiplier tubes. This class is publicly extended by TObservatory.
 */

#ifndef TCamera_hpp
#define TCamera_hpp

#include "TPlane3.hpp"
#include "TSegmentedData.hpp"
#include "THistogramList.hpp"
#include "TRawData.hpp"
#include "TRay.hpp"
#include "TResponseFunction.hpp"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TVirtualFFT.h"
#include <iostream>

class TCamera {
    
private:
    
    Double_t fFocalLength;
    
    Double_t fHeight;
    
    Int_t fNumberTubesY;
    
    Double_t fWidth;
    
    Int_t fNumberTubesX;
    
    Double_t fPMTResolution;
    
    TResponseFunction* fResponseFunction;
    
    Bool_t fCheckBackCollision = false;
    
public:
    
    TCamera(Double_t focalLength, Double_t width, Int_t numberTubesX, Double_t height, Int_t numberTubesY, Double_t PMTResolution, TResponseFunction* responseFunction, Bool_t checkBackCollision);
    
    Double_t FocalLength();
    
    Bool_t CheckCollision(TVector3 position);
    
    TVector3 GetViewDirection(Int_t bin);
    
    TVector2 GetPixelPosition(Int_t bin);
    
    Int_t GetBin(TVector2 pixelPosition);
    
    TSegmentedData SegmentedData(TRawData data);
    
    THistogramList PhotonHistograms(TSegmentedData data);
    
    THistogramList VoltageHistograms(THistogramList photonHistograms, Int_t nFrequencyBins);

};

#endif
