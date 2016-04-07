//
//  TCamera.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TCamera.h"
#include "TVirtualFFT.h"
#include <iostream>

TCamera::TCamera(Double_t focalLength, Double_t width, Int_t numberTubesX, Double_t height, Int_t numberTubesY, Double_t PMTResolution, TResponseFunction* responseFunction, Bool_t checkBackCollision) {
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
    TSegmentedData segmentedData = TSegmentedData(fNumberTubesX * fNumberTubesY);
    for (Int_t i = 0; i < rawData.Size(); i++) {
        if (TMath::Abs(rawData.GetX(i)) > fWidth / 2 || TMath::Abs(rawData.GetY(i)) > fHeight / 2) {
            continue;
        }
        segmentedData.AddPoint(rawData.GetT(i), GetBin(TVector2(rawData.GetX(i), rawData.GetY(i))));
    }
    return segmentedData;
}

THistogramArray TCamera::PhotonHistograms(TSegmentedData parsedData) {
    THistogramArray photonHistograms = THistogramArray();
    Double_t minTime = parsedData.GetMinTime();
    Double_t maxTime = parsedData.GetMaxTime();
    Int_t nHistoBins = (maxTime - minTime) / fPMTResolution;
    for (Int_t bin = 0; bin < parsedData.GetNBins(); bin++) {
        if(parsedData.GetSegment(bin)->size() > 0) {
            TVector2 position = GetPixelPosition(bin);
            Double_t x = position.X();
            Double_t y = position.Y();
            TH1D histogram = TH1D(Form("pmt-x%f-y%f", x, y), Form("Photomultiplier Tube at x = %f, y = %f", x, y), nHistoBins, minTime, maxTime);
            for (Double_t time: *parsedData.GetSegment(bin)) {
                histogram.Fill(time);
            }
            photonHistograms.AddHistogram(x, y, histogram);
        }
    }
    return photonHistograms;
}

THistogramArray TCamera::VoltageHistograms(THistogramArray photonHistograms, Int_t nFrequencyBins) {
    THistogramArray voltageOutput = THistogramArray();
    for (std::list<TPixelData>::iterator iter = photonHistograms.Begin(); iter != photonHistograms.End(); iter++) {
        Double_t test = fResponseFunction->ResponseTime();
        
        std::cout << &fResponseFunction << std::endl;
        
        TH1* tFluxHisto = 0;
        TH1D fluxHisto = *iter;
        tFluxHisto = fluxHisto.FFT(tFluxHisto, "MAG");
        Double_t tFluxReal[nFrequencyBins];
        Double_t tFluxComp[nFrequencyBins];
        
        std::cout << &fResponseFunction << std::endl;
        
        // There is some problem with this line.
        TVirtualFFT::GetCurrentTransform()->GetPointsComplex(tFluxReal, tFluxComp);

        TH1* tResponseHisto = 0;
        TH1D responseHistogram;
        std::cout << &fResponseFunction << std::endl;
        Double_t responseTime = fResponseFunction->ResponseTime();
        responseHistogram = fResponseFunction->ResponseHistogram(fPMTResolution);
        tResponseHisto = responseHistogram.FFT(tResponseHisto, "MAG");
        Double_t tResponseReal[nFrequencyBins];
        Double_t tResponseComp[nFrequencyBins];
        TVirtualFFT::GetCurrentTransform()->GetPointsComplex(tResponseReal, tResponseComp);
        Double_t tProductReal[nFrequencyBins];
        Double_t tProductComp[nFrequencyBins];
        for (Int_t i = 0; i < nFrequencyBins; i++) {
//            tProductReal[i] = tFluxReal[i] * tProductReal[i] - tFluxComp[i] * tProductComp[i];
//            tProductComp[i] = tFluxReal[i] * tProductComp[i] + tFluxComp[i] * tProductReal[i];
        }
        TVirtualFFT *reverseTransform = TVirtualFFT::FFT(1, &nFrequencyBins, "C2R M");
        reverseTransform->SetPointsComplex(tProductReal, tProductComp);
        reverseTransform->Transform();
        TH1D* result = 0;
        TH1::TransformHisto(reverseTransform, result, "Re");
        voltageOutput.AddHistogram((*iter).X(), (*iter).Y(), *result);
        
        delete tFluxHisto;
        delete tResponseHisto;
    }
    return voltageOutput;
}

Int_t TCamera::GetBin(TVector2 position) {
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