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
#include "TAnalysis.h"

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

Int_t TCamera::GetBin(Double_t x, Double_t y) {
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

/*
 * Use the least squares method on the equation derived for time as a function of angle. To find the angle for each data point, find the angle between its position vector and a vector in both the x-y plane and the shower-detector plane. First project it into the the shower-pixel plane.
 */
std::vector<Double_t>* TCamera::ReconstructShower(TSegmentedData data, TTelescope telescope) {
    TPlane3 showerPlane = EstimateShowerPlane(data, telescope);
    TVector3 groundIntersection =  showerPlane.IntersectWithXYPlane();
    std::vector<Double_t> averages = std::vector<Double_t>();
    std::vector<Long_t> weight = std::vector<Long_t>();
    std::vector<Double_t> angles = std::vector<Double_t>();
    for (Int_t i = 0; i < data.GetNBins(); i++) {
        averages.push_back(TAnalysis::SumArray(*data.GetSegment(i)) / data.GetSegment(i)->size());
        weight.push_back(data.GetSegment(i)->size());
        angles.push_back(GetOutwardDirection(telescope, i).Angle(groundIntersection));
    }
    Double_t bestImpactParam = 0;
    Double_t impactRange[] = {0, 100000};
    Double_t impactStepSize = 100;
    Double_t bestShowerAngle = 0;
    Double_t angleStep = 0.1;
    Double_t bestSquare = 1e100;
    for (Double_t impactParam = impactRange[0]; impactParam < impactRange[1]; impactParam += impactStepSize) {
        for (Double_t showerAngle = 0; showerAngle < TMath::Pi(); showerAngle += angleStep) {
            Double_t squareSum = 0;
            Double_t impactAngle = TMath::Pi() / 2 - showerAngle;
            Int_t nearestIndex = 0;
            for (Int_t i = 0; i < angles.size(); i++) {
                if (TMath::Abs(angles[i] - impactAngle) < TMath::Abs(angles[nearestIndex] - impactAngle)) {
                    nearestIndex = i;
                }
            }
            double_t t0 = angles[nearestIndex];
            for (Int_t i = 0; i < averages.size(); i++) {
                squareSum += (t0 + impactParam / TRay::fLightSpeed * TMath::Tan((TMath::Pi() - showerAngle - angles[i]) / 2) - averages[i]) * (t0 + impactParam / TRay::fLightSpeed * TMath::Tan((TMath::Pi() - showerAngle - angles[i]) / 2) - averages[i]);
            }
            if (squareSum < bestSquare) {
                bestImpactParam = impactParam;
                bestShowerAngle = showerAngle;
                bestSquare = squareSum;
            }
        }
    }
    std::vector<Double_t>* output = new std::vector<Double_t>();
    output->push_back(bestImpactParam);
    output->push_back(bestShowerAngle);
    return output;
}

/*
 * Find values of a, b, and c in the plane equation. Propagate each pixel's direction out to the same amount. For given values of a, b, and c, let d = a * xTelescope + b * yTelescope + c * zTelescope. For this plane equation, use the least squares method on the distance between each point and the plane.  Weight squares by the number of photons observed by the pixel.
 */
TPlane3 TCamera::EstimateShowerPlane(TSegmentedData data, TTelescope telescope) {
    std::vector<TVector3> direction = std::vector<TVector3>();
    std::vector<Long_t> weight = std::vector<Long_t>();
    for (Int_t i = 0; i < data.GetNBins(); i++) {
        direction.push_back(GetOutwardDirection(telescope, i));
        weight.push_back(data.GetSegment(i)->size());
    }
    // Find the coefficients of the equation.
    TVector3 centerOfCurvature = telescope.GetCenterOfCurvature();
    TPlane3 bestPlane = TPlane3(TVector3(1, 0, 0), centerOfCurvature);
    Double_t angleStep = 1;
    Double_t bestSquare = 1e100;
    for (Double_t theta = 0; theta <= TMath::Pi(); theta += angleStep) {
        for (Double_t phi = -TMath::Pi(); phi <= TMath::Pi(); phi += angleStep) {
            Double_t a = TMath::Sin(theta) * TMath::Cos(phi);
            Double_t b = TMath::Sin(theta) * TMath::Sin(phi);
            Double_t c = TMath::Cos(theta);
            TPlane3 plane = TPlane3(TVector3(a, b, c), centerOfCurvature);
            Double_t squareSum = 0;
            for (Int_t i = 0; i < direction.size(); i++) {
                squareSum += plane.ShortestDistance(direction[i]) * weight[i];
            }
            if (squareSum < bestSquare) {
                bestSquare = squareSum;
                bestPlane = plane;
            }
        }
    }
    return bestPlane;
}