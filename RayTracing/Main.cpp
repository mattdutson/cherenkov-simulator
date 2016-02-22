/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class is used for testing the functionality of TTelescope, TRay, and TPlane3.
 */

#include "TTelescope.h"
#include "TAnalysis.h"
#include "TFile.h"
#include <iostream>

using namespace std;

void CollectRMSData();

void TestPointImage();

int main(int argc, const char* argv[]) {
    CollectRMSData();
    TestPointImage();
}

void CollectRMSData() {
    TFile file("/Users/Matthew/Documents/XCode/RayTracing/Output/rms-output.root", "RECREATE");
    
    // Sets the dimensions of the TProfile
    Int_t nBinsX = 500;
    Double_t xLow = 0;
    Double_t xUp = TMath::Pi() / 6;
    Double_t yLow = 0;
    Double_t yUp = 0.2;
    
    // Sets the dimensions of the path
    Double_t zDistance = 20000;
    Double_t minAngle = 0;
    Double_t maxAngle = TMath::Pi() / 6;
    
    // Sets the number of data points collected
    Double_t timeDelay = 1e-8;
    Int_t sampleNumber = 100;
    
    // Sets the properties of the mirror
    Double_t radius = 6;
    Short_t mirrorTypes[] = {0, 1};
    Double_t fNumbers[] = {0.5, 1.0, 1.5, 2.0};
    Double_t focalLengths[] = {2.75, 2.8, 2.85, 2.9, 2.95, 3.0};
    
    // Run the simulations
    std::vector<Double_t> RMS = *new std::vector<Double_t>();
    std::vector<Double_t> angle = *new std::vector<Double_t>();
    for (Short_t mirrorType: mirrorTypes) {
        for (Double_t fNumber: fNumbers) {
            for (Double_t focalLength: focalLengths) {
                
                // Generate data points
                TTelescope telescope = *new TTelescope(0, mirrorType, radius, focalLength, fNumber);
                TAnalysis::FindRMSVsAngle(RMS, angle, telescope, sampleNumber, timeDelay, minAngle, maxAngle, zDistance);
                
                // Format the graph title and name
                TString name;
                TString title;
                if (mirrorType == 0) {
                    name = Form("sphr-%f-%f", fNumber, focalLength);
                    title = Form("Spherical Mirror, F-Number: %f, Focal Length: %f", fNumber, focalLength);
                }
                else {
                    name = Form("para-%f-%f", fNumber, focalLength);
                    title = Form("Parabolic Mirror, F-Number: %f, Focal Length: %f", fNumber, focalLength);
                }
                
                // Create the profile and format the axes
                TProfile profile = *new TProfile(name, title, nBinsX, xLow, xUp, yLow, yUp);
                profile.GetXaxis()->SetTitle("Angle (rad)");
                profile.GetYaxis()->SetTitle("RMS Deviation of Distance (m)");
                TAnalysis::FillProfile(angle, RMS, profile);
                profile.Write(name);
            }
        }
    }
    file.Close();
}

void TestPointImage() {
    TFile file("/Users/Matthew/Documents/XCode/RayTracing/Output/point-test.root", "RECREATE");
    
    // Sets the dimensions of the histogram
    Int_t nBinsX = 1000;
    Double_t xLow = -0.1;
    Double_t xUp = 0.1;
    Int_t nBinsY = 1000;
    Double_t yLow = -0.1;
    Double_t yUp = 0.1;
    
    // Sets the points being observed
    Double_t zDistance = 20000;
    Double_t heights[] = {0, 1000, 2000, 3000};
    
    // Sets the number of data points collected
    Int_t sampleNumber = 1000;
    
    //Sets the properties of the mirror
    Short_t mirrorType = 0;
    Double_t radius = 6;
    Double_t focalLength = 3;
    Double_t fNumber = 1;
    
    // Run the simulations
    TTelescope telescope = *new TTelescope(0, mirrorType, radius ,focalLength, fNumber);
    std::vector<Double_t> x = *new std::vector<Double_t>();
    std::vector<Double_t> y = *new std::vector<Double_t>();
    for (Double_t height: heights) {
        
        // Generate data points
        telescope.ViewPoint(*new TVector3(height, 0, zDistance), sampleNumber, x, y);
        
        // Format the graph title and name
        TString name = Form("dist-%f-height-%f", zDistance, height);
        TString title = Form("Image of Point at Distance %f and Height %f", zDistance, height);
        
        // Create a histogram and fill it with data points
        TH2D histogram = *new TH2D(name, title, nBinsX, xLow, xUp, nBinsY, yLow, yUp);
        histogram.GetXaxis()->SetTitle("x (meters)");
        histogram.GetYaxis()->SetTitle("y (meters)");
        TAnalysis::FillHistogram(x, y, histogram);
        histogram.Write(name);
    }
    file.Close();
}