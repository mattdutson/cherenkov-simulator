//
//  TCamera.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/23/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TCamera.h"
#include "TFile.h"
#include "TH1D.h"

Double_t TCamera::GetX(Int_t bin) {
    if (bin >= fNumberTubesX * fNumberTubesY) {
        throw std::invalid_argument("");
    }
    else {
        Int_t xIndex = bin % fNumberTubesX;
        return xIndex * fWidth / (Double_t) fNumberTubesX - fWidth / 2.0;
    }
}

Double_t TCamera::GetY(Int_t bin) {
    if (bin >= fNumberTubesX * fNumberTubesY) {
        throw std::invalid_argument("");
    }
    else {
        Int_t yIndex = (bin - bin % fNumberTubesX) / fNumberTubesX;
        return yIndex * fHeight / (Double_t) fNumberTubesY - fHeight / 2.0;
    }
}

Int_t TCamera::GetBin(Int_t x, Int_t y) {
    Int_t xBin = (x + fWidth / 2) / fWidth * fNumberTubesX;
    Int_t yBin = (y + fHeight / 2) / fHeight * fNumberTubesY;
    return yBin * fNumberTubesX + xBin;
}

TCamera::TCamera() {
}

TCamera::TCamera(Double_t height, Int_t numberTubesY, Double_t width, Int_t numberTubesX, Double_t PMTResponseTime, bool transparent) {
    fTransparent = transparent;
    fHeight = height;
    fNumberTubesY = numberTubesY;
    fWidth = width;
    fNumberTubesX = numberTubesX;
    fPMTResponseTime = PMTResponseTime;
}

TSegmentedData TCamera::ParseData(TRawData rawData) {
    TSegmentedData parsedData = TSegmentedData(fNumberTubesX * fNumberTubesY);
    for (int i = 0; i < rawData.Size(); i++) {
        Double_t x = rawData.GetX(i);
        Double_t y = rawData.GetY(i);
        Double_t t = rawData.GetT(i);
        if (TMath::Abs(x) > fWidth / 2 || TMath::Abs(y) > fHeight / 2) {
            continue;
        }
        else {
            parsedData.AddPoint(t, GetBin(x, y));
        }
    }
    return parsedData;
}

void TCamera::WriteDataToFile(TString filename, TSegmentedData parsedData) {
    Double_t minTime = parsedData.GetMinTime();
    Double_t maxTime = parsedData.GetMaxTime();
    TFile file(filename, "RECREATE");
    for (int bin = 0; bin < fNumberTubesX * fNumberTubesY; bin++) {
        Int_t nBinsx = (maxTime - minTime) / fPMTResponseTime;
        Double_t x = GetX(bin);
        Double_t y = GetY(bin);
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

bool TCamera::CheckCollision(TVector3 position) {
    if (fTransparent) {
        return false;
    }
    else if (TMath::Abs(position.X()) > (fWidth / 2.0) || TMath::Abs(position.Y()) > (fHeight / 2.0)) {
        return false;
    }
    else {
        return true;
    }
}