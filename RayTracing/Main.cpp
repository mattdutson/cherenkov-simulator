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

void testViewShower() {
    TFile file("/Users/Matthew/Documents/XCode/RayTracing/Output/view_shower_output.root", "RECREATE");
    
    TTelescope telescope1 = *new TTelescope(5, 2.5, 0, 1, 1, 50);
    TTelescope telescope2 = *new TTelescope(10, 5, 0, 2, 2, 100);
    TTelescope telescope3 = *new TTelescope(6, 3, 0, 3, 3, 100);
    
    TRay shower1 = *new TRay(*new TVector3(1000, 0, 1000), *new TVector3(0, 0, -1));
    TRay shower2 = *new TRay(*new TVector3(10000, 0, 1000), *new TVector3(0, 1, -1));
    TRay shower3 = *new TRay(*new TVector3(10000, 0, 1000), *new TVector3(0, 0, -1));
    TRay shower4 = *new TRay(*new TVector3(20000, 0, 3000), *new TVector3(0, 0, -1));
    
    TGraph graph1_1 = telescope1.ViewShower(shower1, 1e-9);
    graph1_1.Write("graph1_1");

    TGraph graph1_2 = telescope1.ViewShower(shower2, 1e-9);
    graph1_2.Write("graph1_2");

    TGraph graph1_3 = telescope1.ViewShower(shower3, 1e-9);
    graph1_3.Write("graph1_3");
    
    TGraph graph2_1 = telescope2.ViewShower(shower1, 1e-9);
    graph2_1.Write("graph2_1");
    
    TGraph graph2_2 = telescope2.ViewShower(shower2, 1e-9);
    graph2_2.Write("graph2_2");
    
    TGraph graph2_3 = telescope2.ViewShower(shower3, 1e-9);
    graph2_3.Write("graph2_3");
    
    file.Close();
}

int main(int argc, const char* argv[]) {
    testViewShower();
    return 1;
}
