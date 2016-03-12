//
//  TMirror.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/12/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TMirror.h"

TMirror::TMirror(Short_t mirrorShape, Short_t mirrorType, Double_t radius, Double_t focalLength, Double_t fNumber, Double_t inclination, Double_t azimuth, TVector3 centerOfCurvature) {
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
    
    fRadius = radius;
    fFocalLength = focalLength;
    fCrossDiameter = focalLength / fNumber;
    fCenterOfCurvature = centerOfCurvature;
    
    // Initialize angles
    fInclination = inclination;
    fAzimuth = azimuth;
    
    // Initialize the axis of the mirror
    fMirrorAxis = TVector3(0, 0, 1);
    RotateIn(fMirrorAxis);
}

TVector3* TMirror::GetMirrorImpact() {
    Double_t xRandom = 0;
    Double_t yRandom = 0;
    
    // Select a random point in a circle
    if (fMirrorShape == 0) {
        bool iterate = true;
        while (iterate) {
            xRandom = (fRandom->Rndm() - 0.5) * fCrossDiameter;
            yRandom = (fRandom->Rndm() - 0.5) * fCrossDiameter;
            if (xRandom * xRandom + yRandom * yRandom <= fCrossDiameter * fCrossDiameter / 4) {
                iterate = false;
            }
        }
    }
    
    // Select a random point in a square
    else if (fMirrorShape == 1) {
        xRandom = (fRandom->Rndm() - 0.5) * fCrossDiameter;
        yRandom = (fRandom->Rndm() - 0.5) * fCrossDiameter;
    }
    
    // Find the z-component corresponding to xRandom and yRandom based on the equation of the mirror
    TVector3 relativePosition;
    if (fMirrorType == 0) {
        relativePosition = TVector3(xRandom, yRandom, -fRadius + (fRadius - TMath::Sqrt(fRadius * fRadius - xRandom * xRandom - yRandom * yRandom)));
    }
    else if (fMirrorType == 1) {
        relativePosition = TVector3(xRandom, yRandom, -fRadius + ((xRandom * xRandom / (2 * fRadius) + yRandom * yRandom / (2 * fRadius))));
    }
    
    // Rotate and translate the relative position into the correct frame
    RotateIn(relativePosition);
    return new TVector3(fCenterOfCurvature + relativePosition);
}

TVector3* TMirror::GetMirrorNormal(TVector3 mirrorImpact) {
    TVector3 mirrorNormal;
    
    // For circular mirrors, the normal vector points directly to the center of curvature
    if (fMirrorType == 0) {
        mirrorNormal = fCenterOfCurvature - mirrorImpact;
    }
    
    // For parabolic mirrors, the normal vector is found from the equation z = x^2 / (2R) + y^2 / (2R)
    else if (fMirrorType == 1) {
        mirrorImpact = mirrorImpact - fCenterOfCurvature;
        RotateOut(mirrorImpact);
        mirrorNormal = TVector3(-mirrorImpact.X(), -mirrorImpact.Y(), fRadius);
        RotateIn(mirrorNormal);
    }
    return new TVector3(mirrorNormal);
}

void TMirror::RotateIn(TVector3& vector) {
    vector.RotateY(fInclination);
    vector.RotateX(fAzimuth);
}

void TMirror::RotateOut(TVector3& vector) {
    vector.RotateY(-fInclination);
    vector.RotateX(-fAzimuth);
}

void TMirror::TranslateIn(TVector3& vector) {
    vector += fCenterOfCurvature;
}

void TMirror::TranslateOut(TVector3& vector) {
    vector -= fCenterOfCurvature;
}