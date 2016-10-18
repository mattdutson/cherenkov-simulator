/*
 * CherenkovSimulator - TCoordinates.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TRawData.
 */

#include "TRawData.hpp"

TRawData::TRawData() {
    fXData = std::vector<Double_t>();
    fYData = std::vector<Double_t>();
    fTData = std::vector<Double_t>();
}

void TRawData::PushBack(Double_t x, Double_t y, Double_t t) {
    fXData.push_back(x);
    fYData.push_back(y);
    fTData.push_back(t);
}

void TRawData::Clear() {
    fXData.clear();
    fYData.clear();
    fTData.clear();
}

std::vector<Double_t> TRawData::GetXData() {
    return fXData;
}

std::vector<Double_t> TRawData::GetYData() {
    return fYData;
}

std::vector<Double_t> TRawData::GetTData() {
    return fTData;
}

Double_t TRawData::GetX(Int_t index) {
    return fXData[index];
}

Double_t TRawData::GetY(Int_t index) {
    return fYData[index];
}

Double_t TRawData::GetT(Int_t index) {
    return fTData[index];
}

Double_t TRawData::Size() {
    return fXData.size();
}