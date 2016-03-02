//
//  TDataCollection.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 3/1/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TDataCollection.h"

TDataCollection::TDataCollection() {
    fXData = std::vector<Double_t>();
    fYData = std::vector<Double_t>();
    fTData = std::vector<Double_t>();
    fXMin = 1e100;
    fXMax = -1e100;
    fYMin = 1e100;
    fYMax = -1e100;
    fTMin = 1e100;
    fTMax = -1e100;
}

void TDataCollection::PushBack(Double_t x, Double_t y, Double_t t) {
    fXData.push_back(x);
    fYData.push_back(y);
    fTData.push_back(t);
    if (x < fXMin) {
        fXMin = x;
    }
    if (x > fXMax) {
        fXMax = x;
    }
    if (y < fYMin) {
        fYMin = y;
    }
    if (y > fYMax) {
        fYMax = y;
    }
    if (t < fTMin) {
        fTMin = t;
    }
    if (t > fTMin) {
        fTMax = t;
    }
}

void TDataCollection::Clear() {
    fXData.clear();
    fYData.clear();
    fTData.clear();
}

std::vector<Double_t> TDataCollection::GetXData() {
    return fXData;
}

std::vector<Double_t> TDataCollection::GetYData() {
    return fYData;
}

std::vector<Double_t> TDataCollection::GetTData() {
    return fTData;
}