/*
 * Created by Matthew Dutson on 1/29/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This file contains the implementation of "TTelescope.h". See the header file for method descriptions.
 */

#include "TTelescope.h"
#include <iostream>

TTelescope::TTelescope(Double_t radius, Double_t focalLength, Double_t inclinationAngle, Double_t height, Double_t width, Double_t groudHeight) {
    
    // Initialize miscellaneous variables
    fHeight = height;
    fWidth = width;
    fRadius = radius;
    fInclination = inclinationAngle;

    // Initializes the axis of the mirror
    fMirrorAxis = *new TVector3(1, 0, 0);
    fMirrorAxis.RotateY(inclinationAngle);
    
    // Initializes the focal plane using the inclination angle, the focal length, and the focal plane distance
    TVector3 focalPlaneCenter = *new TVector3(-fRadius + focalLength, 0, 0);
    focalPlaneCenter.RotateY(inclinationAngle);
    fFocalPlane = *new TPlane3(fMirrorAxis, focalPlaneCenter);
    
    // Initializes the ground plane using the height above ground and the assumption that the ground is flat
    fGroundPlane = *new TPlane3(*new TVector3(0, 0, 1), *new TVector3(0, 0, -groudHeight));
}

TVector3 TTelescope::RayDetectionByMirror(TRay shower) {
    
    // Decides on a random point where the ray will strike the mirror
    Double_t yRandom = (fRandom->Rndm() - 0.5) * fWidth;
    Double_t zRandom = (fRandom->Rndm() - 0.5) * fHeight;
    
    // Finds where that point is located on the tangent plane
    TVector3 backPlanePoisition = *new TVector3(-fRadius, yRandom, zRandom);
    backPlanePoisition.RotateY(fInclination);
    
    // Finds the point on the mirror corresponding the back plane position and creates the detected rays
    TVector3 mirrorImpact = backPlanePoisition + fMirrorAxis.Unit() * (fRadius - sqrt(fRadius * fRadius - yRandom * yRandom - zRandom * zRandom));
    TRay detectedRay = *new TRay(mirrorImpact, mirrorImpact - shower.GetPosition());
    
    // Reflects the ray from the mirror and propagates it to the pixel plane
    detectedRay.ReflectFromPlane(*new TPlane3(-mirrorImpact, mirrorImpact));
    detectedRay.PropagateToPlane(fFocalPlane);
    return detectedRay.GetPosition();
}

TGraph TTelescope::ViewShower(TRay shower, Double_t timeDelay) {
    
    // Creates arrays to store the output data
    Int_t numberOfPoints = (Int_t) ((shower.TimeToPlane(fGroundPlane)) / timeDelay) + 2;
    Double_t y[numberOfPoints];
    Double_t z[numberOfPoints];

    // Steps the shower along its path and runs the ray detection algorithm at each point
    for(Int_t i = 0; i < numberOfPoints; i++) {
        TVector3 planeDetection = RayDetectionByMirror(shower);
        planeDetection.RotateY(-fInclination);
        y[i] = planeDetection.Y();
        z[i] = planeDetection.Z();
        shower.IncrementPosition(timeDelay);
    }
    
    // Return the graph we produced
    return *new TGraph(numberOfPoints, y, z);
}