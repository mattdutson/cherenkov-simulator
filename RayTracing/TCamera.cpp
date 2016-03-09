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

/*
 * Use the least squares method on the equation derived for time as a function of angle. To find the angle for each data point, find the angle between its position vector and a vector in both the x-y plane and the shower-detector plane. First project it into the the shower-pixel plane.
 */
std::vector<Double_t>* TCamera::ReconstructShower(TSegmentedData data) {
    Double_t impactParam = 0;
    Double_t planeAngle = 0;
    Double_t showerAngle = 0;
    std::vector<Double_t>* output = new std::vector<Double_t>();
    output->push_back(impactParam);
    output->push_back(planeAngle);
    output->push_back(showerAngle);
    return output;
}

Double_t TCamera::GetPixelAzimuthalAngle(TTelescope telescope, Int_t pixel) {
    return GetOutwardDirection(telescope, pixel).GetPosition().Phi();
}

Double_t TCamera::GetPixelElevationAngle(TTelescope telescope, Int_t pixel) {
    return TMath::Pi() / 2 - GetOutwardDirection(telescope, pixel).GetPosition().Theta();
}

TRay TCamera::GetOutwardDirection(TTelescope telescope, Int_t pixel) {
    TVector3 pixelPosition = TVector3(GetX(pixel), GetY(pixel), -telescope.GetRadius() + telescope.GetFocalLength());
    telescope.RotateIn(pixelPosition);
    telescope.TranslateIn(pixelPosition);
    TRay outwardRay = TRay(0, pixelPosition, (telescope.GetCenterOfCurvature() - telescope.GetAxis().Unit() * telescope.GetRadius()) - pixelPosition);
    outwardRay.ReflectFromPlane(TPlane3(telescope.GetAxis(), TVector3(0, 0, 0)));
    outwardRay.IncrementPosition(1);
    return outwardRay;
}

/*
 * Find values of a, b, and c in the plane equation. Propagate each pixel's direction out to the same amount. For given values of a, b, and c, let d = a * xTelescope + b * yTelescope + c * zTelescope. For this plane equation, use the least squares method on the distance between each point and the plane.  Weight squares by the number of photons observed by the pixel.
 */
TPlane3 TCamera::EstimateShowerPlane(TSegmentedData data, TTelescope telescope) {
    std::vector<Double_t> azimuth = std::vector<Double_t>();
    std::vector<Double_t> elevation = std::vector<Double_t>();
    for (Int_t i = 0; i < data.GetNBins(); i++) {
        if (data.GetSegment(i)->size() > 0) {
            azimuth.push_back(GetPixelAzimuthalAngle(telescope, i));
            elevation.push_back(GetPixelElevationAngle(telescope, i));
        }
    }
    // Find the coefficients of the equation: azimuth = a * elevation + b
    Double_t bestA = 0;
    Double_t aRange[] = {-100, 100};
    Double_t aStep = 0.1;
    Double_t bestB = 0;
    Double_t bRange[] = {-100, 100};
    Double_t bStep = 0.1;
    Double_t bestSquareSum = 1e100;
    for (Double_t a = aRange[0]; a <= aRange[1]; a += aStep) {
        for (Double_t b = bRange[0]; b <= bRange[1]; b += bStep) {
            Double_t squareSum = 0;
            for (Int_t i = 0; i < azimuth.size(); i++) {
                squareSum += (azimuth[i] - elevation[i] - b) * (azimuth[i] - a * elevation[i] - b);
            }
            if (squareSum <= bestSquareSum) {
                bestSquareSum = squareSum;
                bestA = a;
                bestB = b;
            }
        }
    }
    
}