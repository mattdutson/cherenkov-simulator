//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#include "TPlane3.h"

TPlane3::TPlane3() {
    fNormal = *new TVector3();
    fD = 0;
}

TPlane3::TPlane3(TVector3 normal, TVector3 point) {
    fNormal = normal;
    fD = normal.Dot(point);
}

TVector3 TPlane3::GetNormal() {
    return fNormal;
}

Double_t TPlane3::GetEquationCoefficient() {
    return fD;
}
