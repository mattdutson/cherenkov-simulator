/*
 * CherenkovSimulator - TConstantIntensity.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 *
 */

#ifndef TConstantIntensity_h
#define TConstantIntensity_h

#include <stdio.h>
#include "T3DScalarFunction.h"

class TConstantIntensity: public T3DScalarFunction {
    
private:
    
    Int_t fSampleRate;
    
    
public:
    
    TConstantIntensity(Int_t fSampleRate);
    
    Int_t GetIntensity(TVector3 position, Double_t time);
    
};

#endif /* TConstantIntensity_h */
