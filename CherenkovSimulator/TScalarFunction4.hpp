/*
 * CherenkovSimulator - TScalarFunction3.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * An abstract class which defines the behavior of a three-dimensional scalar function. This class may be replaced with a built-in Root class.
 */

#ifndef TScalarFunction4_hpp
#define TScalarFunction4_hpp

#include <stdio.h>
#include "TVector3.h"

class T3DScalarFunction {
    
public:
    
    virtual Int_t GetIntensity(TVector3 position, Double_t time) = 0;
    
};

#endif
