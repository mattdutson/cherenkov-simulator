//
//  TCamera.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TCamera.h"

TCamera::TCamera() {}

TCamera::TCamera(Double_t focalLength, Double_t width, Int_t numberTubesX, Double_t height, Int_t numberTubesY, Double_t PMTResolution, TResponseFunction responseFunction, Bool_t checkBackCollision) {
    fFocalLength = focalLength;
    fWidth = width;
    fNumberTubesX = numberTubesX;
    fHeight = height;
    fNumberTubesY = numberTubesY;
    fCheckBackCollision = checkBackCollision;
    fPMTResolution = PMTResolution;
    fResponseFunction = responseFunction;
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

TSegmentedData TCamera::SegmentedData(TRawData rawData) {
    TSegmentedData parsedData = TSegmentedData(fNumberTubesX * fNumberTubesY);
    for (Int_t i = 0; i < rawData.Size(); i++) {
        if (TMath::Abs(rawData.GetX(i)) > fWidth / 2 || TMath::Abs(rawData.GetY(i)) > fHeight / 2) {
            continue;
        }
        parsedData.AddPoint(rawData.GetT(i), GetBin(TVector3(rawData.GetX(i), rawData.GetY(i), 0)));
    }
    return parsedData;
}

THistogramArray TCamera::PixelHistograms(TSegmentedData parsedData) {
    THistogramArray pixelHistograms = THistogramArray(parsedData.GetNBins());
    Double_t minTime = parsedData.GetMinTime();
    Double_t maxTime = parsedData.GetMaxTime();
    Int_t nHistoBins = (maxTime - minTime) / fPMTResolution;
    for (Int_t bin = 0; bin < parsedData.GetNBins(); bin++) {
        TVector2 position = GetPixelPosition(bin);
        Double_t x = position.X();
        Double_t y = position.Y();
        TH1D histogram = TH1D(Form("pmt-x%f-y%f", x, y), Form("Photomultiplier Tube at x = %f, y = %f", x, y), nHistoBins, minTime, maxTime);
        for (Double_t time: *parsedData.GetSegment(bin)) {
            histogram.Fill(time);
        }
        pixelHistograms.SetHistogram(bin, histogram);
    }
    return pixelHistograms;
}

THistogramArray TCamera::VoltageOutput(THistogramArray histograms, Int_t nFrequencyBins) {
    THistogramArray voltageOutput = THistogramArray(histograms.GetNBins());
    TH1D responseHistogram = fResponseFunction.ResponseHistogram(fPMTResolution);
    //
    
    // Convolve the response function with the photon flux.
    for (Int_t bin = 0; bin < histograms.GetNBins(); bin++) {
        TH1* tFluxReal = new TH1D("", "", nFrequencyBins, - 1 / (2 * fPMTResolution), 1 / (2 * fPMTResolution));
        TH1* tFluxComp = new TH1D("", "", nFrequencyBins, - 1 / (2 * fPMTResolution), 1 / (2 * fPMTResolution));
        TH1* tResponseReal = new TH1D("", "", nFrequencyBins, - 1 / (2 * fPMTResolution), 1 / (2 * fPMTResolution));
        TH1* tResponseComp = new TH1D("", "", nFrequencyBins, - 1 / (2 * fPMTResolution), 1 / (2 * fPMTResolution));
        histograms.GetHistogram(bin).FFT(tFluxReal, "RE R2C M");
        histograms.GetHistogram(bin).FFT(tFluxComp, "IM R2C M");
        responseHistogram.FFT(tResponseReal, "RE R2C M");
        responseHistogram.FFT(tResponseComp, "IM R2C M");
        TH1 tProductReal = TH1D("", "", nFrequencyBins, - 1 / (2 * fPMTResolution), 1 / (2 * fPMTResolution));
        TH1 tProductComp = TH1D("", "", nFrequencyBins, - 1 / (2 * fPMTResolution), 1 / (2 * fPMTResolution));
        
        for (Int_t i = 0; i < nFrequencyBins; i++) {
            tProductReal.SetBinContent(i, tFluxReal->GetBinContent(i) * tResponseReal->GetBinContent(i) - tFluxComp->GetBinContent(i) * tResponseComp->GetBinContent(i));
            tProductComp.SetBinContent(i, tFluxReal->GetBinContent(i) * tResponseComp->GetBinContent(i) + tFluxComp->GetBinContent(i) * tResponseReal->GetBinContent(i));
        }
    }
    return voltageOutput;
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