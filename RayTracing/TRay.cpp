//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#include "TMath.h"
#include "TRay.h"

TRay::TRay(TVector3 position, TVector3 direction) {
    fPosition = position;
    fDirection = direction;
}

void TRay::PropagateToPlane(TPlane3 plane) {
    fPosition += (plane.GetEquationCoefficient() - plane.GetNormal().Dot(fPosition)) / plane.GetNormal().Dot(fDirection) * fPosition;
}

void TRay::ReflectFromPlane(TPlane3 plane) {
    fDirection.Rotate(TMath::Pi(), plane.GetNormal());
    fDirection = -fDirection;
}

TVector3 TRay::GetPosition() {
    return fPosition;
}
