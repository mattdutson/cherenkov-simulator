//
//  TAnalysis.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/8/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TAnalysis.h"
#include "TMath.h"

Double_t TAnalysis::FindRMSFromAverage(std::vector<Double_t> yArray, std::vector<Double_t> zArray) {
    VerifyArraySize(yArray, zArray);
    Long_t n = yArray.size();
    Double_t yAverage = SumArray(yArray) / n;
    Double_t zAverage = SumArray(zArray) / n;
    Double_t variance = 0;
    for (Int_t i = 0; i < n; i++) {
        Double_t distanceSquared = TMath::Power(yArray[i] - yAverage, 2) + TMath::Power(zArray[i] - zAverage, 2);
        variance += distanceSquared;
    }
    variance = variance / n;
    Double_t sigma = TMath::Sqrt(variance);
    return sigma;
}
TH2D TAnalysis::MakeDetectionHistogram(std::vector<Double_t> yArray, std::vector<Double_t> zArray, TString title, Int_t nXBins, Double_t xMin, Double_t xMax, Int_t nYBins, Double_t yMin, Double_t yMax) {
    VerifyArraySize(yArray, zArray);
    Long_t n = yArray.size();
    TH2D histogram = *new TH2D(title, title, nXBins, xMin, xMax, nYBins, yMin, yMax);
    
    for (Int_t i = 0; i < n; i++) {
        histogram.Fill(yArray[i], zArray[i]);
    }
    return histogram;
}

TGraph TAnalysis::MakeGraph(std::vector<Double_t> yArray, std::vector<Double_t> zArray) {
    VerifyArraySize(yArray, zArray);
    TGraph graph = *new TGraph((Int_t) yArray.size(), yArray.data(), zArray.data());
    return graph;
}

TGraph TAnalysis::PlotRMSVAngle(TTelescope telescope, Int_t sampleNumber, Double_t timeDelay, Double_t minAngle, Double_t maxAngle, Double_t xDistance) {
    if (minAngle > maxAngle) {
        Double_t temp = minAngle;
        minAngle = maxAngle;
        maxAngle = temp;
    }
    else if (minAngle == maxAngle) {
        throw std::invalid_argument("The minimum and the maximum angles must be different");
    }
    
    std::vector<Double_t> RMS = *new std::vector<Double_t>();
    std::vector<Double_t> angle = *new std::vector<Double_t>();
    Double_t startingHeight = xDistance * TMath::Tan(minAngle);
    Double_t endingHeight = xDistance * TMath::Tan(maxAngle);
    Int_t nSteps = (Int_t) ((endingHeight - startingHeight) / (timeDelay * 3e8)) + 1;
    
    TRay shower = *new TRay(*new TVector3(xDistance, 0, startingHeight), *new TVector3(0, 0, 1));
    
    std::vector<Double_t> y = *new std::vector<Double_t>();
    std::vector<Double_t> z = *new std::vector<Double_t>();
    
    for (Int_t i = 0; i < nSteps; i++) {
        telescope.ViewPoint(shower.GetPosition(), sampleNumber, y, z);
        Double_t standardDeviation = FindRMSFromAverage(y, z);
        RMS.push_back(standardDeviation);
        angle.push_back(TMath::PiOver2() - shower.GetPosition().Theta());
        shower.IncrementPosition(timeDelay);
        y.clear();
        z.clear();
    }
    
    return MakeGraph(angle, RMS);
}

void TAnalysis::VerifyArraySize(std::vector<Double_t> yArray, std::vector<Double_t> zArray) {
    if (yArray.size() != zArray.size()) {
        throw new std::invalid_argument("Both input arrays must be the same size");
    }
}

Double_t TAnalysis::SumArray(std::vector<Double_t> array) {
    Double_t sum = 0;
    for (Double_t d: array) {
        sum += d;
    }
    return sum;
}