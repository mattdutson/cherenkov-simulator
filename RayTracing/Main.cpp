/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class is used for testing the functionality of TTelescope, TRay, and TPlane3.
 */

#include "TTelescope.h"
#include "TAnalysis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>

using namespace std;

void OutputGraph(TGraph graph, TString title, TString writeName) {
    graph.SetDrawOption("AP");
    graph.SetTitle(title);
    graph.Write(writeName);
}

void testViewShower() {
    TFile file("/Users/Matthew/Documents/XCode/RayTracing/Output/view_shower_output.root", "RECREATE");
    
    TTelescope telescope1 = *new TTelescope(1, 5, 2.5, 0, 1, 50);
    TTelescope telescope2 = *new TTelescope(1, 10, 5, 0, 2, 100);
    
    // Recommended by Prof. Bergman
    TTelescope telescope3 = *new TTelescope(1, 6, 3, 0, 3, 100);
    
    TRay shower1 = *new TRay(*new TVector3(1000, 0, 1000), *new TVector3(0, 0, -1));
    TRay shower2 = *new TRay(*new TVector3(10000, 0, 1000), *new TVector3(0, 1, -1));
    TRay shower3 = *new TRay(*new TVector3(10000, 0, 1000), *new TVector3(0, 0, -1));
    
    // Recommended by Prof. Bergman
    TRay shower4 = *new TRay(*new TVector3(20000, 0, 3000), *new TVector3(0, 0, -1));
    
    // Tests the optical properties of the mirror over the entire shower
    std::vector<Double_t> y = *new std::vector<Double_t>();
    std::vector<Double_t> z = *new std::vector<Double_t>();
    
    telescope1.ViewShower(shower1, 1e-9, 10, y, z);
    TH2D histogram1_1 = TAnalysis::MakeDetectionHistogram(y, z, "histogram1_1", 100, -0.1, 0.1, 100, -1, 1);
    histogram1_1.Write();

    telescope1.ViewShower(shower2, 1e-9, 10, y, z);
    TH2D histogram1_2 = TAnalysis::MakeDetectionHistogram(y, z, "histogram1_2", 100, -0.1, 0.1, 100, -1, 1);
    histogram1_2.Write();

    telescope1.ViewShower(shower3, 1e-9, 10, y, z);
    TH2D histogram1_3 = TAnalysis::MakeDetectionHistogram(y, z, "histogram1_3", 100, -0.1, 0.1, 100, -1, 1);
    histogram1_3.Write();

    telescope2.ViewShower(shower1, 1e-9, 10, y, z);
    TH2D histogram2_1 = TAnalysis::MakeDetectionHistogram(y, z, "histogram2_1", 100, -0.1, 0.1, 100, -1, 1);
    histogram2_1.Write();
    
    telescope2.ViewShower(shower2, 1e-9, 10, y, z);
    TH2D histogram2_2 = TAnalysis::MakeDetectionHistogram(y, z, "histogram2_2", 100, -0.1, 0.1, 100, -1, 1);
    histogram2_2.Write();
    
    telescope2.ViewShower(shower3, 1e-9, 10, y, z);
    TH2D histogram2_3 = TAnalysis::MakeDetectionHistogram(y, z, "histogram2_3", 100, -0.1, 0.1, 100, -1, 1);
    histogram2_3.Write();
    
    // Tests the optical properties of the mirror for a single point
    y.clear();
    z.clear();
    telescope3.ViewPoint(*new TVector3(20000, 0, 20000), 200000, y, z);
    TH2D histogram3_0 = TAnalysis::MakeDetectionHistogram(y, z, "histogram3_0", 200, -2, 2, 200, -10, 0);
    histogram3_0.Write();
    
    y.clear();
    z.clear();
    telescope3.ViewPoint(*new TVector3(20000, 0, 3000), 200000, y, z);
    TH2D histogram3_1 = TAnalysis::MakeDetectionHistogram(y, z, "histogram3_1", 200, -0.25, 0.25, 200, -0.8, -0.3);
    histogram3_1.Write();
    
    y.clear();
    z.clear();
    telescope3.ViewPoint(*new TVector3(20000, 0, 2000), 200000, y, z);
    TH2D histogram3_2 = TAnalysis::MakeDetectionHistogram(y, z, "histogram3_2", 200, -0.1, 0.1, 200, -0.5, -0.25);
    histogram3_2.Write();

    y.clear();
    z.clear();
    telescope3.ViewPoint(*new TVector3(20000, 0, 1000), 200000, y, z);
    TH2D histogram3_3 = TAnalysis::MakeDetectionHistogram(y, z, "histogram3_3", 200, -0.1, 0.1, 200, -0.3, -0.1);
    histogram3_3.Write();
    
    y.clear();
    z.clear();
    telescope3.ViewPoint(*new TVector3(20000, 0, 0), 200000, y, z);
    TH2D histogram3_4 = TAnalysis::MakeDetectionHistogram(y, z, "histogram3_4", 200, -0.1, 0.1, 200, -0.1, 0.1);
    histogram3_4.Write();
    
    TGraph RMSVAngle = TAnalysis::PlotRMSVAngle(telescope3, 500, 1e-8, 0, TMath::Pi() /8, 20000);
    RMSVAngle.Write();
    
    file.Close();
}

int main(int argc, const char* argv[]) {
    testViewShower();
    return 1;
}
