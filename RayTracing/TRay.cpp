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
    fDirection = direction;
}

void TRay::PropagateToPlane(TPlane3 plane) {
    
    // This is t in the parametrization of the ray's path as dx = a*dt, dy = b*dt, dz = c*dt
    Double_t t = plane.GetEquationCoefficient() - plane.GetNormal().Dot(fPosition) / plane.GetNormal().Dot(fDirection);
    
    // Check that the ray actually encounters the plane
    if (t <= 0) {
        throw new invalid_argument("The ray never reaches the plane");
    }
    
    // Update the ray's position
    fPosition += (plane.GetEquationCoefficient() - plane.GetNormal().Dot(fPosition)) / plane.GetNormal().Dot(fDirection) * fDirection;
}

void TRay::ReflectFromPlane(TPlane3 plane) {
    
    // Rotate fDirection by pi around the normal vector and reverse it
    fDirection.Rotate(TMath::Pi(), plane.GetNormal());
    fDirection = -fDirection;
}

TVector3 TRay::GetPosition() {
    return fPosition;
}
