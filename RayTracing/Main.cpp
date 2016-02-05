/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class is used for testing the functionality of TTelescope, TRay, and TPlane3.
 */

#include "TTelescope.h"
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
    
    TTelescope telescope1 = *new TTelescope(5, 2.5, 0, 1, 1, 50);
    TTelescope telescope2 = *new TTelescope(10, 5, 0, 2, 2, 100);
    
    // Recommended by Prof. Bergman
    TTelescope telescope3 = *new TTelescope(6, 3, 0, 3, 3, 100);
    
    TRay shower1 = *new TRay(*new TVector3(1000, 0, 1000), *new TVector3(0, 0, -1));
    TRay shower2 = *new TRay(*new TVector3(10000, 0, 1000), *new TVector3(0, 1, -1));
    TRay shower3 = *new TRay(*new TVector3(10000, 0, 1000), *new TVector3(0, 0, -1));
    
    // Recommended by Prof. Bergman
    TRay shower4 = *new TRay(*new TVector3(20000, 0, 3000), *new TVector3(0, 0, -1));
    
    // Tests the optical properties of the mirror over the entire shower
    TGraph shower1_1 = telescope1.ViewShower(shower1, 1e-9, 10);
    OutputGraph(shower1_1, "Small Telescope with Vertical Nearby Shower", "shower1_1");

    TGraph shower1_2 = telescope1.ViewShower(shower2, 1e-9, 10);
    OutputGraph(shower1_2, "Small Telescope with Angled Faraway Shower", "shower1_2");

    TGraph shower1_3 = telescope1.ViewShower(shower3, 1e-9, 10);
    OutputGraph(shower1_3, "Small Telescope with Vertical Faraway Shower", "shower1_3");
    
    TGraph shower2_1 = telescope2.ViewShower(shower1, 1e-9, 10);
    OutputGraph(shower2_1, "Large Telescope with Vertical Nearby Shower", "shower2_1");
    
    TGraph shower2_2 = telescope2.ViewShower(shower2, 1e-9, 10);
    OutputGraph(shower2_2, "Large Telescope with Angled Faraway Shower", "shower2_2");
    
    TGraph shower2_3 = telescope2.ViewShower(shower3, 1e-9, 10);
    OutputGraph(shower2_3, "Large Telescope with Vertical Faraway Shower", "shower2_3");
    
    // Tests the optical properties of the mirror for a single point
    TGraph point4_1 = telescope3.ViewPoint(*new TVector3(20000, 0, 3000), 200000);
    OutputGraph(point4_1, "Medium Telescope with Point at High Altitude", "point4_1");
    
    TGraph point4_2 = telescope3.ViewPoint(*new TVector3(20000, 0, 2000), 200000);
    OutputGraph(point4_2, "Medium Telescope with Point at Medium Altitude", "point4_2");
    
    TGraph point4_3 = telescope3.ViewPoint(*new TVector3(20000, 0, 1000), 200000);
    OutputGraph(point4_3, "Medium Telescope with Point at Low Altitude", "point4_3");
    
    TGraph point4_4 = telescope3.ViewPoint(*new TVector3(20000, 0, 0), 200000);
    OutputGraph(point4_4, "Medium Telescope with Point at Zero Altitude", "point4_4");
    
    file.Close();
}

int main(int argc, const char* argv[]) {
    testViewShower();
    return 1;
}
