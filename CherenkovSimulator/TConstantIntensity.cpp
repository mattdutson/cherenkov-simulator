/*
 * CherenkovSimulator - TCoordinates.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TConstantIntensity.
 */

#include "TConstantIntensity.hpp"

TConstantIntensity::TConstantIntensity(Int_t sampleRate) {
    fSampleRate = sampleRate;
}

Int_t TConstantIntensity::GetIntensity(TVector3 position, Double_t time) {
    return fSampleRate;
}