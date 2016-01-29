/*
 * Created by Matthew Dutson on 1/29/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This file contains the implementation of "TTelescope.h". See the header file for method descriptions.
 */

#include "TTelescope.h"

TTelescope::TTelescope(Double_t focalLength, Double_t inclination, Double_t height, Double_t width, Double_t focalPlaneDistance, Double_t heightAboveGround) {
    
    // Initializes miscellaneous variables
    fRotationAxis = *new TVector3(1, 0, 0);
    fRandom = new TRandom1();
    fHeight = height;
    fWidth = width;
    fFocalLength = focalLength;
    fRadius = focalLength * 2;
    fInclination = inclination;

    // Initializes the axis of the mirror
    fMirrorAxis = *new TVector3(0, 1, 0);
    fMirrorAxis.Rotate(inclination, fRotationAxis);
    
    // Initializes the pixel plane using the inclination angle, the focal length, and the focal plane distance
    TVector3 pixelPlaneCenter = *new TVector3(0, - fRadius + focalPlaneDistance, 0);
    pixelPlaneCenter.Rotate(inclination, fRotationAxis);
    fPixelPlane = *new TPlane3(fMirrorAxis, pixelPlaneCenter);
    
    // Initializes the mirror tangent plane using the inclination angle and the focal length
    TVector3 mirrorBackCenter = *new TVector3(0, -fRadius, 0);
    mirrorBackCenter.Rotate(inclination, fRotationAxis);
    fMirrorBackPlane = *new TPlane3(fMirrorAxis, mirrorBackCenter);
    
    // Initializes the ground plane using the height above ground and the assumption that the ground is flat
    fGroundPlane = *new TPlane3(*new TVector3(0, 0, 1), *new TVector3(0, 0, -heightAboveGround));
}

TVector3 TTelescope::RayDetectionByMirror(TRay shower) {
    
    // Decides on a random point where the ray will strike the mirror
    Double_t zRandom = (fRandom->Rndm() - 0.5) * fHeight;
    Double_t xRandom = (fRandom->Rndm() - 0.5) * fWidth;
    
    // Finds where that point is located on the tangent plane
    TVector3 backPlanePoisition = *new TVector3(xRandom, 0, zRandom);
    backPlanePoisition.Rotate(fInclination, fRotationAxis);
    
    // Finds the point on the mirror corresponding the back plane position and creates the detected rays
    TVector3 mirrorImpact = backPlanePoisition + fMirrorAxis.Unit() * (fRadius - sqrt(fRadius * fRadius - zRandom * zRandom - xRandom * xRandom));
    TRay detectedRay = *new TRay(mirrorImpact, mirrorImpact - shower.GetPosition());
    
    // Reflects the ray from the mirror and propagates it to the pixel plane
    detectedRay.ReflectFromPlane(*new TPlane3(-mirrorImpact, mirrorImpact));
    detectedRay.PropagateToPlane(fPixelPlane);
    return detectedRay.GetPosition();
}

TGraph TTelescope::ViewShower(TRay shower, Double_t timeDelay) {
    
    // Creates arrays to store the output data
    Int_t numberOfPoints = (int) (shower.TimeToPlane(fGroundPlane) / timeDelay) + 2;
    Double_t x[numberOfPoints];
    Double_t z[numberOfPoints];

    // Steps the shower along its path and runs the ray detection algorithm at each point
    for(Int_t i = 0; shower.TimeToPlane(fGroundPlane) > timeDelay; i++) {
        TVector3 planeDetection = RayDetectionByMirror(shower);
        planeDetection.Rotate(-fInclination, fRotationAxis);
        x[i] = planeDetection.x();
        z[i] = planeDetection.z();
        shower.IncrementPosition(timeDelay);
    }
    
    // When the shower gets close to the ground, propagate it to the ground and run the ray detection algorithm one last time
    shower.PropagateToPlane(fGroundPlane);
    TVector3 planeDetection = RayDetectionByMirror(shower);
    planeDetection.Rotate(-fInclination, fRotationAxis);
    x[numberOfPoints - 1] = planeDetection.x();
    z[numberOfPoints - 1] = planeDetection.z();
    
    // Return the graph we produced
    return *new TGraph(numberOfPoints, x, z);
}