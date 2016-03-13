/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class is used for testing the functionality of TTelescope, TRay, and TPlane3.
 */

#include "TTelescope.h"
#include "TAnalysis.h"
#include "TFile.h"
#include "TConstantIntensity.h"
#include <iostream>

using namespace std;

void CollectRMSData();

void TestPointImage();

void TestCameraFunction();

void TestShowerReconstruction();

int main(int argc, const char* argv[]) {
//    CollectRMSData();
//    TestPointImage();
//    TestCameraFunction();
    TestShowerReconstruction();
}

void CollectRMSData() {
    TFile file("/Users/Matthew/Documents/XCode/CherenkovSimulator/Output/rms-output.root", "RECREATE");
    
    // Sets the dimensions of the TProfile
    Int_t nBinsX = 100;
    Double_t xLow = 0;
    Double_t xUp = TMath::Pi() / 6;
    Double_t yLow = 0;
    Double_t yUp = 0.25;
    
    // Sets the dimensions of the path
    Double_t zDistance = 200000;
    Double_t minAngle = 0;
    Double_t maxAngle = TMath::Pi() / 6;
    
    // Sets the number of data points collected
    Double_t timeDelay = 5e-8;
    Int_t sampleNumber = 300;
    
    // Sets the properties of the mirror
    Double_t radius = 6;
    Short_t mirrorTypes[] = {0, 1};
    Double_t fNumbers[] = {1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
    Double_t focalPercentages[] = {95, 96, 97, 98, 99, 100};
    TCamera camera = TCamera(2, 50, 2, 50, 1e-10, false);
    
    // Run the simulations
    std::vector<Double_t> RMS = std::vector<Double_t>();
    std::vector<Double_t> angle = std::vector<Double_t>();
    for (Short_t mirrorType: mirrorTypes) {
        for (Double_t fNumber: fNumbers) {
            for (Double_t focalPercentage: focalPercentages) {
                
                // Generate data points
                TTelescope telescope(0, mirrorType, radius, radius / 200 * focalPercentage, fNumber, &camera);
                TAnalysis::FindRMSVsAngle(RMS, angle, telescope, sampleNumber, timeDelay, minAngle, maxAngle, zDistance);
                
                // Format the graph title and name
                TString name;
                TString title;
                if (mirrorType == 0) {
                    name = Form("sphr-%f-%f", fNumber, focalPercentage);
                    title = Form("Spherical Mirror, F-Number: %f, Focal Length Percentage: %f", fNumber, focalPercentage);
                }
                else {
                    name = Form("para-%f-%f", fNumber, focalPercentage);
                    title = Form("Parabolic Mirror, F-Number: %f, Focal Length Percentage: %f", fNumber, focalPercentage);
                }
                
                // Create the profile and format the axes
                TProfile profile = TProfile(name, title, nBinsX, xLow, xUp, yLow, yUp);
                profile.GetXaxis()->SetTitle("Angle (rad)");
                profile.GetXaxis()->SetRangeUser(xLow, xUp);
                profile.GetYaxis()->SetTitle("RMS Deviation of Distance (m)");
                profile.GetYaxis()->SetRangeUser(yLow, yUp);
                TAnalysis::FillProfile(angle, RMS, profile);
                profile.Write(name);
            }
        }
    }
    file.Close();
}

void TestPointImage() {
    TFile file("/Users/Matthew/Documents/XCode/CherenkovSimulator/Output/point-images.root", "RECREATE");
    
    // Sets the dimensions of the histogram
    Int_t nBinsX = 200;
    Double_t xLow = -0.5;
    Double_t xUp = 0.5;
    Int_t nBinsY = 200;
    Double_t yLow = -0.8;
    Double_t yUp = 0.2;
    
    // Sets the points being observed
    Double_t zDistance = 20000;
    Double_t heights[] = {0, 1000, 2000, 3000};
    
    // Sets the number of data points collected
    Int_t sampleNumber = 50000;
    
    // Sets the properties of the mirror
    Short_t mirrorType = 0;
    Double_t radius = 6;
    Double_t focalLength = 3;
    Double_t fNumber = 1;
    TCamera camera = TCamera(1.5, 50, 1.5, 50, 1e-10, false);
    
    // Run the simulations
    TTelescope telescope(0, mirrorType, radius, focalLength, fNumber, &camera);
    TRawData data = TRawData();
    for (Double_t height: heights) {
        
        // Generate data points
        TConstantIntensity* intensityFunction = new TConstantIntensity(sampleNumber);
        TShower shower = TShower(TRay(0, TVector3(height, 0, zDistance), TVector3()), intensityFunction);
        telescope.ViewPoint(shower, data);
        
        // Format the graph title and name
        TString name = Form("dist-%f-height-%f", zDistance, height);
        TString title = Form("Image of Point at Distance %f and Height %f", zDistance, height);
        
        // Create a histogram and fill it with data points
        TH2D histogram = TH2D(name, title, nBinsX, xLow, xUp, nBinsY, yLow, yUp);
        histogram.GetXaxis()->SetTitle("x (meters)");
        histogram.GetYaxis()->SetTitle("y (meters)");
        TAnalysis::FillHistogram(data.GetYData(), data.GetXData(), histogram);
        histogram.Write(name);
        delete intensityFunction;
    }
    file.Close();
}

void TestCameraFunction() {
    TFile file("/Users/Matthew/Documents/XCode/CherenkovSimulator/Output/shower-path.root", "RECREATE");
    
    // Set up the camera
    TCamera camera = TCamera(2, 50, 2, 50, 1e-7, false);
    
    // Set the properties of the mirror
    Short_t mirrorType = 0;
    Double_t radius = 6;
    Double_t focalLength = 3;
    Double_t fNumber = 1;
    
    // Set the number of points being observed
    Int_t sampleNumber = 100;
    Double_t delayTime = 1e-8;
    
    // Set up the shower
    TConstantIntensity* intensityFunction = new TConstantIntensity(sampleNumber);
    TShower shower = TShower(TRay(0, TVector3(3000, 0, 20000), TVector3(-1, 0, 0)), intensityFunction);
    
    // Set up the telescope
    TTelescope telescope = TTelescope(0, mirrorType, radius, focalLength, fNumber, &camera);
    
    // Arrays to store data
    TRawData data = TRawData();
    
    telescope.ViewShower(shower, delayTime, data);
    TH2D histogram = TH2D("shower-path", "Shower Path (Height: 3000 m)", 50, -1, 1, 50, -1, 1);
    TAnalysis::FillHistogram(data.GetYData(), data.GetXData(), histogram);
    histogram.Write();
    file.Close();
    TSegmentedData parsedData = telescope.GetCamera()->ParseData(data);
    telescope.GetCamera()->WriteDataToFile("/Users/Matthew/Documents/XCode/CherenkovSimulator/Output/camera-data.root", parsedData);
    delete intensityFunction;
}

void TestShowerReconstruction() {
    
    // Set up the camera
    TCamera camera = TCamera(2, 50, 2, 50, 1e-7, false);
    
    // Set the properties of the mirror
    Short_t mirrorType = 0;
    Double_t radius = 6;
    Double_t focalLength = 3;
    Double_t fNumber = 1;
    
    // Set the number of points being observed
    Int_t sampleNumber = 100;
    Double_t delayTime = 1e-8;
    
    // Set up the shower
    TConstantIntensity* intensityFunction = new TConstantIntensity(sampleNumber);
    TShower shower = TShower(TRay(0, TVector3(3000, 0, 20000), TVector3(-1, 0, 0)), intensityFunction);
    
    // Set up the telescope
    TTelescope telescope = TTelescope(0, mirrorType, radius, focalLength, fNumber, &camera);
    
    // Arrays to store data
    TRawData data = TRawData();
    
    telescope.ViewShower(shower, delayTime, data);
    TSegmentedData parsedData = telescope.GetCamera()->ParseData(data);
    std::vector<Double_t> output = *telescope.GetCamera()->ReconstructShower(parsedData, telescope);
    cout << "Impact Parameter: " << output[0] << endl;
    cout << "Shower Angle: " << output[1] << endl;
}