/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class is used for testing the functionality of TTelescope, TRay, and TPlane3.
 */

#include "TTelescope.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>

using namespace std;

int main(int argc, const char* argv[]) {
    TFile f("test.root", "RECREATE");
    TTelescope telescope = *new TTelescope(1, 0, 2, 2, 1, 50);
    TRay shower = *new TRay(*new TVector3(0, 5000, 5000), *new TVector3(0, 0, -1));
    TGraph graph = telescope.ViewShower(shower, 1e-5);
    graph.Write("graph");
    f.Close();
    return 1;
}