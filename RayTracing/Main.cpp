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

int main(int argc, const char* argv[]) {
    TFile file("/Users/Matthew/Documents/XCode/RayTracing/Output/output.root", "RECREATE");
    
    // Recommended by Prof. Bergman
    TTelescope telescope1 = *new TTelescope(0, 0, 6, 3, 1);
    
    // Recommended by Prof. Bergman
    TRay shower1 = *new TRay(*new TVector3(20000, 0, 3000), *new TVector3(0, 0, -1));
    
    file.Close();
}
