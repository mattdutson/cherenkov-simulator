/*
 * CherenkovSimulator - TShower.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Represents a cosmic ray shower. Identical to a TRay, but with the addition of an intensity function which determines the number of photons emitted by the shower.
 */

#ifndef TShower_hpp
#define TShower_hpp

#include "TRay.hpp"
#include "TScalarFunction4.hpp"

class TShower: public TRay {
    
private:
    T3DScalarFunction* fScalarFunction;
    
public:
    TShower(TRay ray, T3DScalarFunction* scalarFunction);
    
    Int_t GetIntensity();
};

#endif /* TShower_h */
