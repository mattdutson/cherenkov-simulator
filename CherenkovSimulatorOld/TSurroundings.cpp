/*
 * CherenkovSimulator - TSurroundings.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TSurroundings.
 */

#include "TSurroundings.hpp"

TSurroundings::TSurroundings() {}

TSurroundings::TSurroundings(TPlane3 groundPlane, TScalarFunction4* absorptionCoefficients, Double_t spatialIntegrationStep) {
    fGroundPlane = groundPlane;
    fAbsorptionCoefficients = absorptionCoefficients;
    fSpatialIntegrationStep = spatialIntegrationStep;
}

TPlane3 TSurroundings::GroundPlane() {
    return fGroundPlane;
}

Double_t TSurroundings::GetDimmingPercentage(TShower shower, TVector3 endingPoint) {
    TRay ray = TRay(shower.GetTime(), shower.GetPosition(), endingPoint - shower.GetPosition());
    Double_t transmitted = 1;
}