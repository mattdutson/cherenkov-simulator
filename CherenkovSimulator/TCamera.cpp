//
//  TCamera.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TCamera.h"

TPlane3 TCamera::FocalPlane() {
    return fFocalPlane;
}

bool TCamera::CheckCollision(TVector3 position) {
    if (TMath::Abs(position.X()) > (fWidth / 2.0) || TMath::Abs(position.Y()) > (fHeight / 2.0)) {
        return false;
    }
    else {
        return true;
    }
}

TSegmentedData TCamera::ParseData(TRawData rawData) {
    TSegmentedData parsedData = TSegmentedData(fNumberTubesX * fNumberTubesY);
    for (Int_t i = 0; i < rawData.Size(); i++) {
        parsedData.AddPoint(rawData.GetT(i), GetBin(TVector3(rawData.GetX(i), rawData.GetY(i), 0)));
    }
    return parsedData;
}

Int_t TCamera::GetBin(TVector3 position) {
    Int_t xBin = (position.X() + fWidth / 2) / fWidth * fNumberTubesX;
    Int_t yBin = (position.Y() + fHeight / 2) / fHeight * fNumberTubesY;
    return yBin * fNumberTubesX + xBin;
}

TVector3 TCamera::GetViewDirection(Int_t bin) {
    Double_t x;
    Double_t y;
    if (bin >= fNumberTubesX * fNumberTubesY) {
        throw std::invalid_argument("");
    }
    else {
        Int_t xBin = bin % fNumberTubesX;
        Int_t yBin = (bin - bin % fNumberTubesX) / fNumberTubesX;
        y = yBin * fHeight / (Double_t) fNumberTubesY - fHeight / 2.0;
        x = xBin * fWidth / (Double_t) fNumberTubesX - fWidth / 2.0;
    }
    return TVector3(-x, -y, fFocalLength).Unit();
}