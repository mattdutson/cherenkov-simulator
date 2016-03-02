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
}

void TDataCollection::PushBack(Double_t x, Double_t y, Double_t t) {
    fXData.push_back(x);
    fYData.push_back(y);
    fTData.push_back(t);
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

Double_t TDataCollection::GetX(Int_t index) {
    return fXData[index];
}

Double_t TDataCollection::GetY(Int_t index) {
    return fYData[index];
}

Double_t TDataCollection::GetT(Int_t index) {
    return fTData[index];
}

Double_t TDataCollection::Size() {
    return fXData.size();
}