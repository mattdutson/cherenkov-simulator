//
//  TSegmentedData.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/8/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TSegmentedData.h"

TSegmentedData::TSegmentedData(Int_t nBins) {
    fSegmentedData = new std::vector<Double_t>*[nBins];
    for (Int_t i = 0; i < nBins; i++) {
        fSegmentedData[i] = new std::vector<Double_t>();
    }
    fNBins = nBins;
    fMinTime = 1e100;
    fMaxTime = 1e-100;
}

void TSegmentedData::AddPoint(Double_t time, Int_t bin) {
    fSegmentedData[bin]->push_back(time);
    if (time > fMaxTime) {
        fMaxTime = time;
    }
    if (time < fMinTime) {
        fMinTime = time;
    }
}

Double_t TSegmentedData::GetMinTime() {
    return fMinTime;
}

Double_t TSegmentedData::GetMaxTime() {
    return fMaxTime;
}

std::vector<Double_t>* TSegmentedData::GetSegment(Int_t bin) {
    return fSegmentedData[bin];
}

Int_t TSegmentedData::GetNBins() {
    return fNBins;
}