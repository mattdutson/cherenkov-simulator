/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class is used for testing the functionality of TTelescope, TRay, and TPlane3.
 */

#include "TShower.h"

TShower::TShower(TRay ray, T3DScalarFunction* scalarFunction): TRay(ray) {
    fScalarFunction = scalarFunction;
}

Int_t TShower::GetIntensity() {
    return fScalarFunction->GetIntensity(fPosition, fTime);
}