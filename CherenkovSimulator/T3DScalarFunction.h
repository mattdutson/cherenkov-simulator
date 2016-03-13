//
//  T3DScalarFunction.hpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/28/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef T3DScalarFunction_h
#define T3DScalarFunction_h

#include <stdio.h>
#include "TVector3.h"

class T3DScalarFunction {
    
public:
    
    virtual Int_t GetIntensity(TVector3 position, Double_t time) = 0;
    
};

#endif /* T3DScalarFunction_h */
