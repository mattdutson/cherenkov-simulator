/*
 * CherenkovSimulator - TCoordinates.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TShower.
 */

#include "TShower.hpp"

TShower::TShower(TRay ray, T3DScalarFunction* scalarFunction): TRay(ray) {
    fScalarFunction = scalarFunction;
}

Int_t TShower::GetIntensity() {
    return fScalarFunction->GetIntensity(fPosition, fTime);
}