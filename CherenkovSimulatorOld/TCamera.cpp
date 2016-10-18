/*
 * CherenkovSimulator - TCamera.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TCamera.
 */

#include "TCamera.hpp"

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

THistogramList TCamera::PhotonHistograms(TSegmentedData parsedData) {
    THistogramList photonHistograms = THistogramList();
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

THistogramList TCamera::VoltageHistograms(THistogramList photonHistograms, Int_t nFrequencyBins) {
    THistogramList voltageOutput = THistogramList();
    
    for (std::list<TPixelData>::iterator iter = photonHistograms.Begin(); iter != photonHistograms.End(); iter++) {
        TH1D photonHistogram = *iter;
        TH1D responseHistogram = fResponseFunction->ResponseHistogram(fPMTResolution);
        
        Int_t nSamplePoints = photonHistogram.GetNbinsX();
        Int_t nResponsePoints = responseHistogram.GetNbinsX();
        
        // This should be modified later to make the number of points some power of 2
        Int_t nPointsMin = nSamplePoints + nResponsePoints;
        Int_t nextPowerOfTwo = ((Int_t) TMath::Log2(nPointsMin)) + 1;
        Int_t powerOfTwo = TMath::Power(2, nextPowerOfTwo);
        Double_t photonData[powerOfTwo];
        Double_t responseData[powerOfTwo];
        
        // Create the arrays and pad them with zeros.
        for(Int_t i = 0; i < powerOfTwo; i++) {
            if (i < nSamplePoints) {
                photonData[i] = photonHistogram.GetBinContent(i);
            }
            else {
                photonData[i] = 0;
            }
            if (i < nResponsePoints) {
                responseData[i] = responseHistogram.GetBinContent(i);
            }
            else {
                responseData[i] = 0;
            }
        }
        
        Int_t nTransformPoints = powerOfTwo + 1;
        
        // Create an FFT object and transform photon data
        TVirtualFFT* photonFFT = TVirtualFFT::FFT(1, &powerOfTwo, "R2C M K");
        photonFFT->SetPoints(photonData);
        photonFFT->Transform();
        Double_t photonFrequenciesRe[nTransformPoints];
        Double_t photonFrequenciesIm[nTransformPoints];
        photonFFT->GetPointsComplex(photonFrequenciesRe, photonFrequenciesIm);

        // Create an FFT object and transform response data
        TVirtualFFT* responseFFT = TVirtualFFT::FFT(1, &powerOfTwo, "R2C M K");
        responseFFT->SetPoints(responseData);
        responseFFT->Transform();
        Double_t responseFrequenciesRe[nTransformPoints];
        Double_t responseFrequenciesIm[nTransformPoints];
        responseFFT->GetPointsComplex(responseFrequenciesRe, responseFrequenciesIm);
        
        // Create the arrays for the product
        Double_t productRe[nTransformPoints];
        Double_t productIm[nTransformPoints];
        
        // Compute the product of photon and response frequencies
        for (Int_t i = 0; i < nTransformPoints; i++) {
            productRe[i] = photonFrequenciesRe[i] * responseFrequenciesRe[i] - photonFrequenciesIm[i] * responseFrequenciesIm[i];
            productIm[i] = photonFrequenciesRe[i] * responseFrequenciesIm[i] + photonFrequenciesIm[i] * responseFrequenciesRe[i];
        }
        
        // Take the reverse FFT
        TVirtualFFT *reverseTransform = TVirtualFFT::FFT(1, &nTransformPoints, "C2R M K");
        reverseTransform->SetPointsComplex(productRe, productIm);
        reverseTransform->Transform();
        Double_t outputSignal[powerOfTwo];
        reverseTransform->GetPoints(outputSignal);
        
        // Make a histogram from the signal data
        TH1D result = TH1D(Form("Signal at x=%f, y=%f", iter->X(), iter->Y()), "Signal", nPointsMin, photonHistogram.GetXaxis()->GetXmin(), photonHistogram.GetXaxis()->GetXmin() +fPMTResolution * nPointsMin);
        for (Int_t i = 0; i < nPointsMin; i++) {
            result.SetBinContent(i, outputSignal[i]);
        }
        
        // Add the resulting histogram to the output
        voltageOutput.AddHistogram((*iter).X(), (*iter).Y(), result);
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