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
    
    fRandom = new TRandom1();
}


TVector3 TTelescope::rayDetectionByMirror(TRay shower) {
    
    Double_t heightMultiplier = fRandom->Rndm() - 0.5;
    
    Double_t widthMultiplier = fRandom->Rndm() - 0.5;
    
    Double_t z = heightMultiplier * getHeight();
    
    Double_t x = widthMultiplier * getWidth();
    
    TVector3 backPlanePoisition = *new TVector3(x, 0, z);
    
    TVector3 rotationAxis = *new TVector3(1, 0, 0);
    
    backPlanePoisition.Rotate(getInclination(), rotationAxis);
    
    TVector3 rayDirection = backPlanePoisition - shower.GetPosition();
    
    TVector3 mirrorAxis = getMirrorAxis();
    Double_t mirrorRadius = getRadius();
    
    TVector3 mirrorImpact = backPlanePoisition + mirrorAxis.Unit() * (mirrorRadius - sqrt(mirrorRadius * mirrorRadius - z * z - x * x));
    
    TRay detectedRay = *new TRay(mirrorImpact, rayDirection);
    
    detectedRay.ReflectFromPlane(*new TPlane3(-mirrorImpact, mirrorImpact));
    
    detectedRay.PropagateToPlane(getFocalPlane());
    
    return detectedRay.GetPosition();
}



TGraph TTelescope::viewShower(TRay shower, Double_t timeDelay) {
    Int_t numberOfPoints = (int) (shower.distanceToPlane(getGroundPlane()) / (3e8) / timeDelay) + 2;
    Double_t x[numberOfPoints];
    Double_t z[numberOfPoints];
    Double_t distanceStep = (3e8) * timeDelay;
    TVector3 rotationAxis = *new TVector3(1, 0, 0);
    Int_t index = 0;
    while (shower.distanceToPlane(getGroundPlane()) > distanceStep) {
        TVector3 planeDetection = rayDetectionByMirror(shower);
        planeDetection.Rotate(-getInclination(), rotationAxis);
        x[index] = planeDetection.X();
        z[index] = planeDetection.Z();
        shower.incrementPosition(timeDelay);
    }
    shower.PropagateToPlane(getGroundPlane());
    TVector3 planeDetection = rayDetectionByMirror(shower);
    planeDetection.Rotate(-getInclination(), rotationAxis);
    x[numberOfPoints - 1] = planeDetection.X();
    z[numberOfPoints - 1] = planeDetection.Z();
    return *new TGraph(numberOfPoints, x, z);
}

Double_t TTelescope::getHeight() {
    return height;
}

Double_t TTelescope::getWidth() {
    return width;
}

Double_t TTelescope::getInclination() {
    return inclination;
}

TVector3 TTelescope::getMirrorAxis() {
    return focalPlane.GetNormal();
}

Double_t TTelescope::getRadius() {
    return focalLength*2;
}

TPlane3 TTelescope::getFocalPlane() {
    return focalPlane;
}

TPlane3 TTelescope::getGroundPlane() {
    return groundPlane;
}