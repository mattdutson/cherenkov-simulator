/*
 * CherenkovSimulator - TCoordinates.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TUtility.
 */

#include "TUtility.hpp"

Double_t TUtility::SumArray(std::vector<Double_t> array) {
    Double_t sum = 0;
    for (Double_t d: array) {
        sum += d;
    }
    return sum;
}

void TUtility::WriteHistogramFile(TString filename, THistogramList histograms) {
    TFile file(filename, "RECREATE");
    for (std::list<TPixelData>::iterator iter = histograms.Begin(); iter != histograms.End(); iter++) {
        (*iter).Write();
    }
    file.Close();
}

void TUtility::FillHistogram(std::vector<Double_t> array1, std::vector<Double_t> array2, TH2D& histogram) {
    Long_t n = array1.size();
    for (Int_t i = 0; i < n; i++) {
        histogram.Fill(array1[i], array2[i]);
    }
}

void TUtility::FillProfile(std::vector<Double_t> array1, std::vector<Double_t> array2, TProfile& profile) {
    Long_t n = array1.size();
    for (Int_t i = 0; i < n; i++) {
        profile.Fill(array1[i], array2[i]);
    }
}