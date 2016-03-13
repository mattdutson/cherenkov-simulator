/*
 * Created by Matthew Dutson on 1/29/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This file contains the implementation of "TTelescope.h". See the header file for method descriptions. A left-handed coordinate system is used, with the z-axis pointing along the axis of the telescope and the x-axis oriented parallel to the horizontal.
 */

#include "TTelescope.h"

void TTelescope::ViewPointPrivate(TShower shower, TRawData& data) {
    
    // Steps the shower along its path and runs the ray detection algorithm at each point
    for(Int_t i = 0; i < shower.GetIntensity(); i++) {
        TRay* planeDetection = RayDetection(shower);
        
        // If the detected ray hit the camera, skip this iteration
        if(planeDetection == nullptr) {
            delete planeDetection;
            continue;
        }
        
        TVector3 position =  planeDetection->GetPosition();
        
        // Change coordinates to the telescope frame and store data in the array
        RotateOut(position);
        TranslateOut(position);
        data.PushBack(position.X(), position.Y(), shower.GetTime());
        delete planeDetection;
    }
}

TRay* TTelescope::RayDetection(TShower shower) {

    // Find where the detected ray hits the mirror and the normal vector at that point
    TVector3* mirrorImpact = GetMirrorImpact();
    TVector3* mirrorNormal = GetMirrorNormal(*mirrorImpact);
    
    // Create the detected ray
    TRay* detectedRay = new TRay(shower.GetTime(), shower.GetPosition(), (*mirrorImpact) - shower.GetPosition());
    
    // Propagate the detected ray to the focal plane and check whether it collided with the mirror
    detectedRay->PropagateToPlane(fFocalPlane);
    TVector3 position = detectedRay->GetPosition();
    RotateOut(position);
    TranslateOut(position);
    
    if(fCamera->CheckCollision(position)) {
        delete mirrorImpact;
        delete mirrorNormal;
        delete detectedRay;
        return nullptr;
    }
    
    // If the detected ray didn't hit the camera, continue with the simulation
    detectedRay->PropagateToPoint(*mirrorImpact);
    
    // Reflects the ray from the mirror and propagates it to the pixel plane
    detectedRay->ReflectFromPlane(TPlane3(*mirrorNormal, *mirrorImpact));
    detectedRay->PropagateToPlane(fFocalPlane);
    delete mirrorImpact;
    delete mirrorNormal;
    return detectedRay;
}

TCamera* TTelescope::GetCamera() {
    return fCamera;
}

TTelescope::TTelescope(TCamera* camera, TMirror mirror): TTelescope(TPlane3(TVector3(1, 0, 0), TVector3(0, 0, 0)), camera, mirror) {}

TTelescope::TTelescope(TPlane3 groundPlane, TCamera* camera, TMirror mirror): TMirror(mirror) {

    // Initialize member variables
    fCamera = camera;
    fGroundPlane = groundPlane;
    
    // Initialize the focal plane
    TVector3 relativePosition(0, 0, -fRadius + fFocalLength);
    RotateIn(relativePosition);
    TVector3 focalPlaneCenter = fCenterOfCurvature + relativePosition;
    fFocalPlane = TPlane3(fMirrorAxis, focalPlaneCenter);
}

void TTelescope::ViewShower(TShower shower, Double_t timeDelay, TRawData& data) {
    
    // Clear the array before starting.
    data.Clear();
    
    // Creates arrays to store the output data
    Int_t numberOfSteps = (Int_t) (((shower.TimeToPlane(fGroundPlane)) / timeDelay) + 2);
    
    // Steps the shower along its path and runs the ray detection algorithm at each point
    for(Int_t i = 0; i < numberOfSteps; i++) {
        ViewPointPrivate(shower, data);
        shower.IncrementPosition(timeDelay);
    }
}

void TTelescope::ViewPoint(TShower shower, TRawData& data) {
    data.Clear();
    ViewPointPrivate(shower, data);
}

TVector3 TTelescope::GetOutwardDirection(Double_t pixelX, Double_t pixelY) {
        TVector3 pixelPosition = TVector3(pixelX, pixelY, -fFocalLength);
        RotateIn(pixelPosition);
        TranslateIn(pixelPosition);
        TRay outwardRay = TRay(0, pixelPosition, (fCenterOfCurvature - fMirrorAxis.Unit() * fRadius) - pixelPosition);
        outwardRay.ReflectFromPlane(TPlane3(fMirrorAxis, TVector3(0, 0, 0)));
        outwardRay.IncrementPosition(1);
        return outwardRay.GetPosition();
}