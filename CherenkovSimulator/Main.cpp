/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class is used for testing the functionality of TTelescope, TRay, and TPlane3.
 */

#include "TObservatory.h"
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
    
    // Run the simulations
    std::vector<Double_t> RMS = std::vector<Double_t>();
    std::vector<Double_t> angle = std::vector<Double_t>();
    for (Short_t mirrorType: mirrorTypes) {
        for (Double_t fNumber: fNumbers) {
            for (Double_t focalPercentage: focalPercentages) {
                
                TMirror mirror = TMirror(mirrorType, 0, radius, radius / 2 / fNumber);
                TCoordinates coordinates = TCoordinates(0, 0, TVector3(0, 0, 0));
                TCamera camera = TCamera(radius / 200 * focalPercentage, 5, 100, 5, 100, 1e-8, true);
                TSurroundings surroundings = TSurroundings(TPlane3(TVector3(0, 1, 0), TVector3(0, 0, 0)));
                
                // Generate data points
                TObservatory observatory(mirror, camera, coordinates, surroundings);
                TAnalysis::FindRMSVsAngle(RMS, angle, observatory, sampleNumber, timeDelay, minAngle, maxAngle, zDistance);
                
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
    
    // Run the simulations
    TMirror mirror = TMirror(mirrorType, 0, radius, radius / 2 / fNumber);
    TCoordinates coordinates = TCoordinates(0, 0, TVector3(0, 0, 0));
    TCamera camera = TCamera(focalLength, 5, 100, 5, 100, 1e-8, false);
    TSurroundings surroundings = TSurroundings(TPlane3(TVector3(0, 1, 0), TVector3(0, 0, 0)));
    
    TObservatory observatory(mirror, camera, coordinates, surroundings);
    for (Double_t height: heights) {
        
        // Generate data points
        TConstantIntensity* intensityFunction = new TConstantIntensity(sampleNumber);
        TShower shower = TShower(TRay(0, TVector3(0, height, zDistance), TVector3(0, -1, 0)), intensityFunction);
        TRawData data = observatory.ViewPoint(shower);
        
        // Format the graph title and name
        TString name = Form("dist-%f-height-%f", zDistance, height);
        TString title = Form("Image of Point at Distance %f and Height %f", zDistance, height);
        
        // Create a histogram and fill it with data points
        TH2D histogram = TH2D(name, title, nBinsX, xLow, xUp, nBinsY, yLow, yUp);
        histogram.GetXaxis()->SetTitle("x (meters)");
        histogram.GetYaxis()->SetTitle("y (meters)");
        TAnalysis::FillHistogram(data.GetXData(), data.GetYData(), histogram);
        histogram.Write(name);
        delete intensityFunction;
    }
    file.Close();
}

void TestCameraFunction() {
    TFile file("/Users/Matthew/Documents/XCode/CherenkovSimulator/Output/shower-path.root", "RECREATE");
    
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
    TShower shower = TShower(TRay(0, TVector3(0, 3000, 20000), TVector3(0, -1, 0)), intensityFunction);
    
    // Set up the telescope
    TMirror mirror = TMirror(mirrorType, 0, radius, radius / 2 / fNumber);
    TCoordinates coordinates = TCoordinates(0, 0, TVector3(0, 0, 0));
    TCamera camera = TCamera(focalLength, 5, 100, 5, 100, 1e-8, false);
    TSurroundings surroundings = TSurroundings(TPlane3(TVector3(0, 1, 0), TVector3(0, 0, 0)));
    TObservatory observatory = TObservatory(mirror, camera, coordinates, surroundings);
    
    TRawData data = observatory.ViewShower(shower, delayTime);
    TH2D histogram = TH2D("shower-path", "Shower Path (Height: 3000 m)", 50, -1, 1, 50, -1, 1);
    TAnalysis::FillHistogram(data.GetXData(), data.GetYData(), histogram);
    histogram.Write();
    file.Close();
    
    // Parse the data and write it to a file.
    TSegmentedData parsedData = observatory.ParseData(data);
    observatory.WriteDataToFile("/Users/Matthew/Documents/XCode/CherenkovSimulator/Output/camera-data.root", parsedData);
    delete intensityFunction;
}

void TestShowerReconstruction() {
    
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
    TShower shower = TShower(TRay(0, TVector3(0, 3000, 20000), TVector3(0, -1, 0)), intensityFunction);
    
    // Set up the observatory
    TMirror mirror = TMirror(mirrorType, 0, radius, radius / 2 / fNumber);
    TCoordinates coordinates = TCoordinates(0, 0, TVector3(0, 0, 0));
    TCamera camera = TCamera(focalLength, 5, 100, 5, 100, 1e-8, false);
    TSurroundings surroundings = TSurroundings(TPlane3(TVector3(0, 1, 0), TVector3(0, 0, 0)));
    TObservatory observatory = TObservatory(mirror, camera, coordinates, surroundings);
    
    TRawData data1 = observatory.ViewShower(shower, delayTime);
    TSegmentedData parsedData1 = observatory.ParseData(data1);
    TRay output = observatory.ReconstructShower(parsedData1);
    std::cout << "Velocity:" << endl;
    output.GetVelocity().Print();
    std::cout << endl;
    std::cout << "Closest Position: " << endl;
    output.GetPosition().Print();
    std::cout << endl;
    
    shower = TShower(TRay(0, TVector3(0, 3000, 20000), TVector3(1, -1, 0)), intensityFunction);
    TRawData data2 = observatory.ViewShower(shower, delayTime);
    TSegmentedData parsedData2 = observatory.ParseData(data2);
    output = observatory.ReconstructShower(parsedData2);
    std::cout << "Velocity:" << endl;
    output.GetVelocity().Print();
    std::cout << endl;
    std::cout << "Closest Position: " << endl;
    output.GetPosition().Print();
}