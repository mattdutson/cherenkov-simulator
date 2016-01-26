//
//  TTelescope.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 1/25/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TTelescope.h"

TTelescope::TTelescope(Double_t focalLength, Double_t inclinationAngle, Double_t height, Double_t width, Double_t focalPlaneDistance, Double_t heightAboveGround) {
    
    groundPlane = *new TPlane3(*new TVector3(0, 0, 1), *new TVector3(0, 0, -heightAboveGround));
    
    TVector3 rotationAxis = *new TVector3(1, 0, 0);
    
    TVector3 mirrorAxis = *new TVector3(0, 1, 0);
    mirrorAxis.Rotate(inclinationAngle, rotationAxis);
    
    TVector3 focalPlaneCenter = *new TVector3(0, - 2 * focalLength + focalPlaneDistance, 0);
    focalPlaneCenter.Rotate(inclinationAngle, rotationAxis);
    
    TVector3 mirrorBackCenter = *new TVector3(0, -2 * focalLength, 0);
    mirrorBackCenter.Rotate(inclinationAngle, rotationAxis);
    
    focalPlane = *new TPlane3(mirrorAxis, focalPlaneCenter);
    mirrorBackPlane = *new TPlane3(mirrorAxis, mirrorBackCenter);
}