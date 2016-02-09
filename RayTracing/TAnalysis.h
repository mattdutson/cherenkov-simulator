//
//  TAnalysis.hpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/8/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TAnalysis_h
#define TAnalysis_h

#include "TH2.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TArrayD.h"
#include "TTelescope.h"
#include <stdio.h>

class TAnalysis {
    
private:
    
    void VerifyArraySize(TArrayD yArray, TArrayD zArray);
    
public:

    TH2D MakeDetectionHistogram(TArrayD yArray, TArrayD zArray, TString title, Int_t nXBins, Double_t xMin, Double_t xMax, Int_t nYBins, Double_t yMin, Double_t yMax);

    TGraph MakeGraph(TArrayD yArray, TArrayD zArray);

    TGraph PlotRMSVAngle(TTelescope telescope, Double_t delayTime, Int_t sampleNumber, Double_t timeDelay, Double_t minAngle, Double_t maxAngle, Double_t xDistance);

    Double_t FindRMSFromAverage(TArrayD yArray, TArrayD zArray);

};

#endif /* TAnalysis_h */
