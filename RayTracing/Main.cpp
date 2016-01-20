/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class is used for testing the functionality of TTelescope, TRay, and TPlane3.
 */

#include <iostream>
#include <math.h>
#include "TTelescope.h"

using namespace std;

int main(int argc, const char *argv[]) {
    
    TTelescope *telescope = new TTelescope(*new TVector3 (0, 0, 0), *new TPlane3(*new TVector3(0, 0, 1), *new TVector3(0, 0, -10)));
    
    // Simulate a ray moving straight up toward the telescope (focal length 1)
    TVector3 source = telescope->RayTrace(*new TVector3(0, 0, 1), *new TVector3(0, 0, 2));
    cout << "Ray moving vertically:" << endl;
    cout << "(" << source.X() << ", " << source.Y() << ", " << source.Z() << ")" << endl;
    
    // Simulate a ray starting at ~45 degree angle (focal length 1)
    source = telescope->RayTrace(*new TVector3(-1, 0, 0.1), *new TVector3(-sqrt(2)/2, 0, -sqrt(2)/2));
    cout << "Ray starting at ~45 degree angle:" << endl;
    cout << "(" << source.X() << ", " << source.Y() << ", " << source.Z() << ")" << endl;    return 0;
    
    return 0;
}
