/*
 * CherenkovSimulator - TConstantIntensity.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Represents a four-dimensional scalar function which is constant over all space and time.
 */

#ifndef TConstantIntensity_hpp
#define TConstantIntensity_hpp

#include "TScalarFunction4.hpp"

#include <stdio.h>

class TConstantIntensity: public TScalarFunction4 {
    
private:
    
    Int_t fSampleRate;
    
    
public:
    
    TConstantIntensity(Int_t fSampleRate);
    
    Int_t GetIntensity(TVector3 position, Double_t time);
    
};

#endif /* TConstantIntensity_h */
