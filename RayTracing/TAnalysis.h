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
    
    static void VerifyArraySize(std::vector<Double_t> yArray, std::vector<Double_t> zArray);
    
    static Double_t SumArray(std::vector<Double_t>);
    
public:

    static TH2D MakeDetectionHistogram(std::vector<Double_t> yArray, std::vector<Double_t> zArray, TString title, Int_t nXBins, Double_t xMin, Double_t xMax, Int_t nYBins, Double_t yMin, Double_t yMax);

    static TGraph MakeGraph(std::vector<Double_t> yArray, std::vector<Double_t> zArray);

    static TGraph PlotRMSVAngle(TTelescope telescope, Int_t sampleNumber, Double_t timeDelay, Double_t minAngle, Double_t maxAngle, Double_t xDistance);

    static Double_t FindRMSFromAverage(std::vector<Double_t> yArray, std::vector<Double_t> zArray);

};

#endif /* TAnalysis_h */
