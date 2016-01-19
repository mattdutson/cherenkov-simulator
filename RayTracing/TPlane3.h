//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#ifndef TPlane3_h
#define TPlane3_h

#include "TVector3.h"

class TPlane3 {
    
private:
    TVector3 fNormal;
    Double_t fD;
    
public:
    TPlane3();
    TPlane3(TVector3 normal, TVector3 point);
    TVector3 GetNormal();
    Double_t GetEquationCoefficient();
};

#endif
