/*
 * Created by Matthew Dutson on 1/29/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This file contains the implementation of "TTelescope.h". See the header file for method descriptions.
 */

#include "TTelescope.h"

void TTelescope::ViewPointPrivate(TVector3 objectPosition, Int_t sampleNumber, std::vector<Double_t>& xArray, std::vector<Double_t>& yArray) {
    
    // Steps the shower along its path and runs the ray detection algorithm at each point
    for(Int_t i = 0; i < sampleNumber; i++) {
        TVector3 planeDetection = RayDetection(objectPosition);
        
        // Change coordinates to the telescope frame and store data in the array
        RotateOut(planeDetection);
        TranslateOut(planeDetection);
        xArray.push_back(planeDetection.X());
        yArray.push_back(planeDetection.Y());
    }
}

TVector3 TTelescope::RayDetection(TVector3 objectPosition) {

    // Find where the detected ray hits the mirror and the normal vector at that point
    TVector3 mirrorImpact = GetMirrorImpact();
    TVector3 mirrorNormal = GetMirrorNormal(mirrorImpact);
    
    // Create the detected ray
    TRay detectedRay = *new TRay(mirrorImpact, mirrorImpact - objectPosition);
    
    // Reflects the ray from the mirror and propagates it to the pixel plane
    detectedRay.ReflectFromPlane(*new TPlane3(mirrorNormal, mirrorImpact));
    detectedRay.PropagateToPlane(fFocalPlane);
    return detectedRay.GetPosition();
}

TVector3 TTelescope::GetMirrorImpact() {
    Double_t xRandom = 0;
    Double_t yRandom = 0;
    
    // Select a random point in a square
    if (fMirrorShape == 0) {
        xRandom = (fRandom->Rndm() - 0.5) * fCrossDiameter;
        yRandom = (fRandom->Rndm() - 0.5) * fCrossDiameter;
    }
    
    // Select a random point in a circle
    else if (fMirrorShape == 1) {
        bool iterate = true;
        while (iterate) {
            xRandom = (fRandom->Rndm() - 0.5) * fCrossDiameter;
            yRandom = (fRandom->Rndm() - 0.5) * fCrossDiameter;
            if (xRandom * xRandom + yRandom * yRandom <= fCrossDiameter * fCrossDiameter / 4) {
                iterate = false;
            }
        }
    }
    
    // Find the z-component corresponding to xRandom and yRandom based on the equation of the mirror
    TVector3 relativePosition;
    if (fMirrorType == 0) {
        relativePosition = *new TVector3(0, 0, -fRadius + (fRadius - TMath::Sqrt(fRadius * fRadius - xRandom * xRandom - yRandom * yRandom)));
    }
    else if (fMirrorType == 1) {
        relativePosition = *new TVector3(0, 0, -fRadius + ((xRandom * xRandom / (2 * fRadius) + yRandom * yRandom / (2 * fRadius))));
    }
    
    // Rotate and translate the relative position into the correct frame
    RotateIn(relativePosition);
    return fCenterOfCurvature + relativePosition;
}

TVector3 TTelescope::GetMirrorNormal(TVector3 mirrorImpact) {
    TVector3 mirrorNormal;
    
    // For circular mirrors, the normal vector points directly to the center of curvature
    if (fMirrorType == 0) {
        mirrorNormal = fCenterOfCurvature - mirrorImpact;
    }

    // For parabolic mirrors, the normal vector is found from the equation z = x^2 / (2R) + y^2 / (2R)
    else if (fMirrorType == 1) {
        mirrorImpact = mirrorImpact - fCenterOfCurvature;
        RotateOut(mirrorImpact);
        mirrorNormal = *new TVector3(-mirrorImpact.X(), -mirrorImpact.Y(), fRadius);
        RotateIn(mirrorNormal);
    }
    return mirrorNormal;
}

void TTelescope::RotateIn(TVector3& vector) {
    vector.RotateY(fInclination);
    vector.RotateX(fAzimuth);
}

void TTelescope::RotateOut(TVector3& vector) {
    vector.RotateY(-fInclination);
    vector.RotateX(-fAzimuth);
}

void TTelescope::TranslateIn(TVector3& vector) {
    vector += fCenterOfCurvature;
}

void TTelescope::TranslateOut(TVector3& vector) {
    vector -= fCenterOfCurvature;
}

TTelescope::TTelescope(Short_t mirrorShape, Short_t mirrorType, Double_t radius, Double_t focalLength, Double_t fNumber): TTelescope(mirrorShape, mirrorType, radius, focalLength, fNumber, 0, 0, *new TVector3(0, 0, 0), *new TPlane3(*new TVector3(1, 0, 0), *new TVector3(0, 0, 0))) {}

TTelescope::TTelescope(Short_t mirrorShape, Short_t mirrorType, Double_t radius, Double_t focalLength, Double_t fNumber, Double_t inclination, Double_t azimuth, TVector3 centerOfCurvature, TPlane3 groundPlane) {
    // Set the mirror shape, checking for invalid input
    if (mirrorShape < 0 || mirrorShape > 1) {
        throw new std::invalid_argument("The mirror shape must lie in the range [0, 1]");
    }
    else {
        fMirrorShape = mirrorShape;
    }
    
    // Set the mirror type, checking for invalid input
    if (mirrorType < 0 || mirrorType > 1) {
        throw new std::invalid_argument("The mirror type must lie in the range [0, 1]");
    }
    else {
        fMirrorType = mirrorType;
    }
    
    // Initialize member variables
    fRadius = radius;
    fCrossDiameter = focalLength / fNumber;
    fCenterOfCurvature = centerOfCurvature;
    fGroundPlane = groundPlane;
    
    // Initialize angles
    fInclination = inclination;
    fAzimuth = azimuth;
    
    // Initialize the axis of the mirror
    fMirrorAxis = *new TVector3(0, 0, 1);
    RotateIn(fMirrorAxis);
    
    // Initialize the focal plane
    TVector3 relativePosition = *new TVector3(0, 0, -fRadius + focalLength);
    RotateIn(relativePosition);
    TVector3 focalPlaneCenter = fCenterOfCurvature + relativePosition;
    fFocalPlane = *new TPlane3(fMirrorAxis, focalPlaneCenter);
}

void TTelescope::ViewShower(TRay shower, Double_t timeDelay, Int_t sampleNumber, std::vector<Double_t>& xArray, std::vector<Double_t>& yArray) {
    
    // Clear the arrays before starting.
    xArray.clear();
    yArray.clear();
    
    // Creates arrays to store the output data
    Int_t numberOfSteps = (Int_t) (((shower.TimeToPlane(fGroundPlane)) / timeDelay) + 2);
    
    // Steps the shower along its path and runs the ray detection algorithm at each point
    for(Int_t i = 0; i < numberOfSteps; i++) {
        ViewPointPrivate(shower.GetPosition(), sampleNumber, xArray, yArray);
        shower.IncrementPosition(timeDelay);
    }
}

void TTelescope::ViewPoint(TVector3 position, Int_t sampleNumber, std::vector<Double_t> &xArray, std::vector<Double_t> &yArray) {
    xArray.clear();
    yArray.clear();
    ViewPointPrivate(position, sampleNumber, xArray, yArray);
}