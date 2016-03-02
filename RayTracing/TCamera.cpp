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

Double_t TCamera::GetPMTX(Int_t xIndex) {
    if (xIndex >= fNumberTubesX) {
        throw std::invalid_argument("");
    }
    else {
        return xIndex * fWidth / (Double_t) fNumberTubesX - fWidth / 2.0;
    }
}

Double_t TCamera::GetPMTY(Int_t yIndex) {
    if (yIndex >= fNumberTubesY) {
        throw std::invalid_argument("");
    }
    else {
        return yIndex * fHeight / (Double_t) fNumberTubesY - fHeight / 2.0;
    }
}

Double_t TCamera::GetXBin(Double_t x) {
    return (x + fWidth / 2) / fWidth * fNumberTubesX;
}

Double_t TCamera::GetYBin(Double_t y) {
    return (y + fHeight / 2) / fHeight * fNumberTubesY;
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
    fMinTime = 1e100;
    fMaxTime = -1e100;
}

std::vector<Double_t>*** TCamera::ParseData(std::vector<Double_t> xArray, std::vector<Double_t> yArray, std::vector<Double_t> timeArray) {
    fMinTime = 1e100;
    fMaxTime = -1e100;
    std::vector<Double_t>*** data = new std::vector<Double_t>**[fNumberTubesX];
    for (int i = 0; i < fNumberTubesX; i++) {
        data[i] = new std::vector<Double_t>*[fNumberTubesY];
        for (int j = 0; j < fNumberTubesY; j++) {
            data[i][j] = new std::vector<Double_t>();
        }
    }
    Long_t n = xArray.size();
    for (int i = 0; i < n; i++) {
        Int_t xBin = GetXBin(xArray[i]);
        Int_t yBin = GetYBin(yArray[i]);
        if (xBin >= fNumberTubesX || yBin >= fNumberTubesY) {
            continue;
        }
        else {
            data[xBin][yBin]->push_back(timeArray[i]);
            if (timeArray[i] > fMaxTime) {
                fMaxTime = timeArray[i];
            }
            if (timeArray[i] < fMinTime) {
                fMinTime = timeArray[i];
            }
        }
    }
    return data;
}

void TCamera::WriteDataToFile(TString filename, std::vector<Double_t> ***data) {
    TFile file(filename, "RECREATE");
    for (int i = 0; i < fNumberTubesX; i++) {
        for (int j = 0; j < fNumberTubesY; j++) {
            Int_t nBinsx = (fMaxTime - fMinTime) / fPMTResponseTime;
            TH1D histogram = TH1D(Form("pmt-x%i-y%i", i, j), Form("Photomultiplier Tube at x = %f, y = %f", GetPMTX(i), GetPMTY(j)), nBinsx, fMinTime, fMaxTime);
            for (Double_t time: *data[i][j]) {
                histogram.Fill(time);
            }
            histogram.Write();
        }
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