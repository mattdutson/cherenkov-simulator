/*
 * Created by Matthew Dutson on 1/29/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This file contains the implementation of "TTelescope.h". See the header file for method descriptions.
 */

#include "TTelescope.h"
#include <iostream>
#include "TMath.h"

TTelescope::TTelescope(Int_t mirrorType, Double_t radius, Double_t focalLength, Double_t inclinationAngle, Double_t size, Double_t groudHeight) {
    
    // Sets the mirror type
    if (mirrorType < 0 || mirrorType > 1) {
        throw new std::invalid_argument("The mirror type must lie in the range [0, 1]");
    }
    else {
        fMirrorType = mirrorType;
    }
    
    // Initialize miscellaneous variables
    fSize = size;
    fRadius = radius;
    fInclination = inclinationAngle;

    // Initializes the axis of the mirror
    fMirrorAxis = *new TVector3(1, 0, 0);
    fMirrorAxis.RotateY(inclinationAngle);
    fMirrorAxis = fMirrorAxis.Unit();
    
    // Initializes the focal plane using the inclination angle, the focal length, and the focal plane distance
    TVector3 focalPlaneCenter = *new TVector3(-fRadius + focalLength, 0, 0);
    focalPlaneCenter.RotateY(inclinationAngle);
    fFocalPlane = *new TPlane3(fMirrorAxis, focalPlaneCenter);
    
    // Initializes the ground plane using the height above ground and the assumption that the ground is flat
    fGroundPlane = *new TPlane3(*new TVector3(0, 0, 1), *new TVector3(0, 0, -groudHeight));
}

TVector3 TTelescope::RayDetection(TVector3 objectPosition) {

    // Finds the point on the mirror corresponding the back plane position and creates the detected rays
    TVector3 mirrorImpact = GetImpactPoint();
    TRay detectedRay = *new TRay(mirrorImpact, mirrorImpact - objectPosition);
    
    // Reflects the ray from the mirror and propagates it to the pixel plane
    detectedRay.ReflectFromPlane(*new TPlane3(-mirrorImpact, mirrorImpact));
    detectedRay.PropagateToPlane(fFocalPlane);
    return detectedRay.GetPosition();
}

void TTelescope::ViewShower(TRay shower, Double_t timeDelay, Int_t sampleNumber, std::vector<Double_t>& yArray, std::vector<Double_t>& zArray) {
    
    if (yArray.size() > 0) {
        yArray.clear();
        zArray.clear();
    }
    
    // Creates arrays to store the output data
    Int_t numberOfSteps = (Int_t) (((shower.TimeToPlane(fGroundPlane)) / timeDelay) + 2);

    // Steps the shower along its path and runs the ray detection algorithm at each point
    for(Int_t i = 0; i < numberOfSteps; i++) {
        ViewPoint(shower.GetPosition(), sampleNumber, yArray, zArray);
        shower.IncrementPosition(timeDelay);
    }
    
}

void TTelescope::ViewPoint(TVector3 objectPosition, Int_t sampleNumber, std::vector<Double_t>& yArray, std::vector<Double_t>& zArray) {
    
    // Steps the shower along its path and runs the ray detection algorithm at each point
    for(Int_t i = 0; i < sampleNumber; i++) {
        TVector3 planeDetection = RayDetection(objectPosition);
        planeDetection.RotateY(-fInclination);
        yArray.push_back(planeDetection.Y());
        zArray.push_back(planeDetection.Z());
    }
}

TVector3 TTelescope::GetImpactPoint() {
    
    TVector3 impactPosition;
    Double_t yRandom;
    Double_t zRandom;
    
    if (fMirrorType == 0) {
        // Decides on a random point where the ray will strike the mirror
        yRandom = (fRandom->Rndm() - 0.5) * fSize;
        zRandom = (fRandom->Rndm() - 0.5) * fSize;
    }
    else if (fMirrorType == 1) {
        bool iterate = true;
        while (iterate) {
            yRandom = (fRandom->Rndm() - 0.5) * fSize;
            zRandom = (fRandom->Rndm() - 0.5) * fSize;
            if (yRandom * yRandom + zRandom * zRandom <= fSize * fSize / 4) {
                iterate = false;
            }
        }
    }
    
    // Finds where that point is located on the tangent plane
    TVector3 backPlanePoisition = *new TVector3(-fRadius, yRandom, zRandom);
    backPlanePoisition.RotateY(fInclination);
    
    // Finds the point on the mirror in front of the random point on the tangent plane
    impactPosition = backPlanePoisition + fMirrorAxis * (fRadius - TMath::Sqrt(fRadius * fRadius - yRandom * yRandom - zRandom * zRandom));
    
    return impactPosition;
}


