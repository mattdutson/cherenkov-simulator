/*
 * CherenkovSimulator - TCoordinates.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TPlane3.
 */

#include "TPlane3.hpp"

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

TVector3 TPlane3::IntersectWithXZPlane() {
    return TVector3((fCoefficient - fNormal.Z()) / fNormal.X(), 0, 1).Unit();
}

TVector3 TPlane3::ProjectOntoPlane(TVector3 point) {
    return point + ShortestDistance(point) * fNormal;
}