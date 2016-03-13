//
//  TShower.hpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/28/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TShower_h
#define TShower_h

#include <stdio.h>
#include "TRay.h"
#include "T3DScalarFunction.h"

class TShower: public TRay {
    
private:
    T3DScalarFunction* fScalarFunction;
    
public:
    TShower(TRay ray, T3DScalarFunction* scalarFunction);
    
    Int_t GetIntensity();
};

#endif /* TShower_h */
