//
//  TCamera.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TCamera.h"

TCamera::TCamera() {}

TCamera::TCamera(Double_t focalLength, Double_t width, Int_t numberTubesX, Double_t height, Int_t numberTubesY, Double_t PMTResolution, Bool_t checkBackCollision) {
    fFocalLength = focalLength;
    fWidth = width;
    fNumberTubesX = numberTubesX;
    fHeight = height;
    fNumberTubesY = numberTubesY;
    fCheckBackCollision = checkBackCollision;
    fPMTResolution = PMTResolution;
}

Double_t TCamera::FocalLength() {
    return fFocalLength;
}

bool TCamera::CheckCollision(TVector3 position) {
    if (!fCheckBackCollision) {
        return false;
    }
    else if (TMath::Abs(position.X()) > (fWidth / 2.0) || TMath::Abs(position.Y()) > (fHeight / 2.0)) {
        return false;
    }
    else {
        return true;
    }
}

TSegmentedData TCamera::ParseData(TRawData rawData) {
    TSegmentedData parsedData = TSegmentedData(fNumberTubesX * fNumberTubesY);
    for (Int_t i = 0; i < rawData.Size(); i++) {
        if (TMath::Abs(rawData.GetX(i)) > fWidth / 2 || TMath::Abs(rawData.GetY(i)) > fHeight / 2) {
            continue;
        }
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
    TVector2 position;
    if (bin >= fNumberTubesX * fNumberTubesY) {
        throw std::invalid_argument("");
    }
    else {
        TVector2 position = GetPixelPosition(bin);
        return TVector3(-position.X(), -position.Y(), fFocalLength).Unit();
    }
    return TVector3();
}

TVector2 TCamera::GetPixelPosition(Int_t bin) {
    Int_t xBin = bin % fNumberTubesX;
    Int_t yBin = (bin - bin % fNumberTubesX) / fNumberTubesX;
    Double_t y = yBin * fHeight / (Double_t) fNumberTubesY - fHeight / 2.0;
    Double_t x = xBin * fWidth / (Double_t) fNumberTubesX - fWidth / 2.0;
    return TVector2(x, y);
}

void TCamera::WriteDataToFile(TString filename, TSegmentedData parsedData) {
    Double_t minTime = parsedData.GetMinTime();
    Double_t maxTime = parsedData.GetMaxTime();
    TFile file(filename, "RECREATE");
    for (int bin = 0; bin < fNumberTubesX * fNumberTubesY; bin++) {
        Int_t nBinsx = (maxTime - minTime) / fPMTResolution;
        TVector2 position = GetPixelPosition(bin);
        Double_t x = position.X();
        Double_t y = position.Y();
        TH1D histogram = TH1D(Form("pmt-x%f-y%f", x, y), Form("Photomultiplier Tube at x = %f, y = %f", x, y), nBinsx, minTime, maxTime);
        if (parsedData.GetSegment(bin)->size() == 0) {
            continue;
        }
        for (Double_t time: *parsedData.GetSegment(bin)) {
            histogram.Fill(time);
        }
        histogram.Write();
    }
    file.Close();
}