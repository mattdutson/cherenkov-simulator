/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This file contains the implementation of "TPlane3.h". See the header file for method descriptions.
 */

#include "TPlane3.h"

TPlane3::TPlane3() {
    fNormal = TVector3();
    fCoefficient = 0;
}

TPlane3::TPlane3(TVector3 normal, TVector3 point) {
    fNormal = normal.Unit();
    
    // The coefficient is calculated by plugging the point into the plane equation
    fCoefficient = normal.Dot(point);
}

TVector3 TPlane3::GetNormal() {
    return fNormal;
}

Double_t TPlane3::GetEquationCoefficient() {
    return fCoefficient;
}

Double_t TPlane3::ShortestDistance(TVector3 point) {
    return (fNormal.Dot(point) - fCoefficient);
}

TVector3 TPlane3::IntersectWithXYPlane() {
    return TVector3((fCoefficient - fNormal.Y()) / fNormal.X(), 1, 0).Unit();
}

TVector3 TPlane3::ProjectOntoPlane(TVector3 point) {
    return point + ShortestDistance(point) * fNormal;
}