/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This file contains the implementation of "TRay.h". See the header file for method descriptions.
 */

#include "TMath.h"
#include "TRay.h"

using namespace std;

TRay::TRay(TVector3 position, TVector3 direction) {
    fPosition = position;
    fVelocity = fLightSpeed * direction.Unit();
}

TVector3 TRay::GetPosition() {
    return fPosition;
}

TVector3 TRay::GetVelocity() {
    return fVelocity;
}

Double_t TRay::TimeToPlane(TPlane3 plane) {
    TVector3 normal = plane.GetNormal();
    Double_t coefficient = plane.GetEquationCoefficient();
    
    // The time it would take for the ray to reach the plane
    Double_t time = (coefficient - normal.Dot(fPosition)) / normal.Dot(fVelocity);
    
    // Check that the ray actually encounters the plane
    if (time < 0) {
        throw new invalid_argument("The ray never reaches the plane");
    }
    
    // If the ray reaches the plane, return the time
    return time;
}

void TRay::PropagateToPlane(TPlane3 plane) {
    IncrementPosition(TimeToPlane(plane));
}

void TRay::ReflectFromPlane(TPlane3 plane) {
    fVelocity -= 2 * fVelocity.Dot(plane.GetNormal()) * plane.GetNormal();
}

void TRay::IncrementPosition(Double_t time) {
    fPosition += fVelocity * time;
}
