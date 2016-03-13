//
//  TShower.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/28/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TShower.h"

TShower::TShower(TRay ray, T3DScalarFunction* scalarFunction): TRay(ray) {
    fScalarFunction = scalarFunction;
}

Int_t TShower::GetIntensity() {
    return fScalarFunction->GetIntensity(fPosition, fTime);
}