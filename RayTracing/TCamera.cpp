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

TCamera::TCamera() {
}

TCamera::TCamera(Double_t height, Int_t numberTubesY, Double_t width, Int_t numberTubesX, Double_t PMTResponseTime) {
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
        Double_t x = xArray[i] + fWidth / 2;
        Double_t y = yArray[i] + fHeight / 2;
        Int_t xBin = (Int_t) (x / fWidth) * fNumberTubesX;
        Int_t yBin = (Int_t) (y / fHeight) * fNumberTubesY;
        data[xBin][yBin]->push_back(timeArray[i]);
        if (timeArray[i] > fMaxTime) {
            fMaxTime = timeArray[i];
        }
        if (timeArray[i] < fMinTime) {
            fMinTime = timeArray[i];
        }
    }
    return data;
}

void TCamera::WriteDataToFile(TString filename, std::vector<Double_t> ***data) {
    TFile file(filename, "RECREATE");
    for (int i = 0; i < fNumberTubesX; i++) {
        for (int j = 0; j < fNumberTubesY; j++) {
            TH1D histogram = TH1D(Form("pmt-x%i-y%i", i, j), Form("Photomultiplier Tube at x = %f, y = %f", i * fWidth / (Double_t) fNumberTubesX - fWidth / 2, j * fHeight / (Double_t) fNumberTubesY - fHeight / 2), (fMaxTime - fMinTime) / fPMTResponseTime, fMinTime, fMaxTime);
            for (Double_t time: *data[i][j]) {
                histogram.Fill(time);
            }
            histogram.Write();
        }
    }
    file.Close();
}