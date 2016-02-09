//
//  TAnalysis.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/8/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TAnalysis.h"
#include "TMath.h"

Double_t TAnalysis::FindRMSFromAverage(TArrayD yArray, TArrayD zArray) {
    VerifyArraySize(yArray, zArray);
    Int_t n = yArray.GetSize();
    Double_t yAverage = yArray.GetSum() / n;
    Double_t zAverage = zArray.GetSum() / n;
    Double_t variance = 0;
    for (Int_t i = 0; i < n; i++) {
        Double_t distanceSquared = TMath::Power(yArray.GetAt(i) - yAverage, 2) + TMath::Power(zArray.GetAt(i) - zAverage, 2);
        variance += distanceSquared;
    }
    variance = variance / n;
    Double_t sigma = TMath::Sqrt(variance);
    return sigma;
}

TH2D TAnalysis::MakeDetectionHistogram(TArrayD yArray, TArrayD zArray, TString title, Int_t nXBins, Double_t xMin, Double_t xMax, Int_t nYBins, Double_t yMin, Double_t yMax) {
    VerifyArraySize(yArray, zArray);
    Int_t n = yArray.GetSize();
    TH2D histogram = *new TH2D(title, title, 100, -1, 1, 100, -1, 1);
    
    for (Int_t i = 0; i < n; i++) {
        histogram.Fill(yArray.GetAt(i), zArray.GetAt(i));
    }
    return histogram;
}

TGraph TAnalysis::MakeGraph(TArrayD yArray, TArrayD zArray) {
    VerifyArraySize(yArray, zArray);
    TGraph graph = *new TGraph(yArray.GetSize(), yArray.GetArray(), zArray.GetArray());
    return graph;
}

TGraph TAnalysis::PlotRMSVAngle(TTelescope telescope, Double_t delayTime, Int_t sampleNumber, Double_t timeDelay, Double_t minAngle, Double_t maxAngle, Double_t xDistance) {
    if (minAngle > maxAngle) {
        Double_t temp = minAngle;
        minAngle = maxAngle;
        maxAngle = temp;
    }
    else if (minAngle == maxAngle) {
        throw std::invalid_argument("The minimum and the maximum angles must be different");
    }
    
    TArrayD RMS = *new TArrayD();
    TArrayD angle = *new TArrayD();
    Double_t startingHeight = xDistance * TMath::Tan(minAngle);
    Double_t endingHeight = xDistance * TMath::Tan(maxAngle);
    Int_t nSteps = (Int_t) ((endingHeight - startingHeight) / timeDelay) + 1;
    
    TRay shower = *new TRay(*new TVector3(xDistance, 0, startingHeight), *new TVector3(0, 0, 1));
    
    TArrayD y = *new TArrayD();
    TArrayD z = *new TArrayD();
    
    for (Int_t i = 0; i < nSteps; i++) {
        telescope.ViewPoint(shower.GetPosition(), sampleNumber, y, z);
        Double_t standardDeviation = FindRMSFromAverage(y, z);
        RMS.SetAt(standardDeviation, i);
        angle.SetAt(TMath::PiOver2() - shower.GetPosition().Theta(), i);
        shower.IncrementPosition(timeDelay);
        y.Reset();
        z.Reset();
    }
    
    return MakeGraph(angle, RMS);
}

void TAnalysis::VerifyArraySize(TArrayD yArray, TArrayD zArray) {
    if (yArray.GetSize() != zArray.GetSize()) {
        throw new std::invalid_argument("Both input arrays must be the same size");
    }
}