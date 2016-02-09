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
    TArrayD y = *new TArrayD();
    TArrayD z = *new TArrayD();
    
    telescope1.ViewShower(shower1, 1e-9, 10, y, z);

    y.Reset();
    z.Reset();
    telescope1.ViewShower(shower2, 1e-9, 10, y, z);

    y.Reset();
    z.Reset();
    telescope1.ViewShower(shower3, 1e-9, 10, y, z);
    
    y.Reset();
    z.Reset();
    telescope2.ViewShower(shower1, 1e-9, 10, y, z);
    
    y.Reset();
    z.Reset();
    telescope2.ViewShower(shower2, 1e-9, 10, y, z);
    
    y.Reset();
    z.Reset();
    telescope2.ViewShower(shower3, 1e-9, 10, y, z);
    
    // Tests the optical properties of the mirror for a single point
    y.Reset();
    z.Reset();
    telescope3.ViewPoint(*new TVector3(20000, 0, 3000), 200000, y, z);
    
    y.Reset();
    z.Reset();
    telescope3.ViewPoint(*new TVector3(20000, 0, 2000), 200000, y, z);
    
    y.Reset();
    z.Reset();
    telescope3.ViewPoint(*new TVector3(20000, 0, 1000), 200000, y, z);
    
    y.Reset();
    z.Reset();
    telescope3.ViewPoint(*new TVector3(20000, 0, 0), 200000, y, z);
    
    file.Close();
}

int main(int argc, const char* argv[]) {
    testViewShower();
    return 1;
}
