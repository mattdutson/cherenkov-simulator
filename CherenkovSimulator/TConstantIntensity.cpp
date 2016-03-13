//
//  TConstantIntensity.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 3/1/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TConstantIntensity.h"
#include "T3DScalarFunction.h"

TConstantIntensity::TConstantIntensity(Int_t sampleRate) {
    fSampleRate = sampleRate;
}

Int_t TConstantIntensity::GetIntensity(TVector3 position, Double_t time) {
    return fSampleRate;
}