/*
 * CherenkovSimulator - TMirror.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TMirror.
 */

#include "TMirror.hpp"
#include "TMath.h"

TMirror::TMirror() {}

TMirror::TMirror(Int_t type, Int_t shape, Double_t radius, Double_t size) {
    fType = type;
    fShape = shape;
    fRadius = radius;
    fSize = size;
}

Double_t TMirror::Radius() {
    return fRadius;
}

TVector3 TMirror::GetMirrorImpact() {
    Double_t xRandom = 0;
    Double_t yRandom = 0;
    if (fShape == 0) {
        Bool_t iterate = true;
        while (iterate) {
            xRandom = (fRng->Rndm() - 0.5) * fSize;
            yRandom = (fRng->Rndm() - 0.5) * fSize;
            if (xRandom * xRandom + yRandom * yRandom <= fSize * fSize / 4) {
                iterate = false;
            }
        }
    }
    else if (fShape == 1) {
        xRandom = (fRng->Rndm() - 0.5) * fSize;
        yRandom = (fRng->Rndm() - 0.5) * fSize;
    }
    TVector3 impactPosition;
    if (fType == 0) {
        impactPosition = TVector3(xRandom, yRandom, -fRadius + (fRadius - TMath::Sqrt(fRadius * fRadius - xRandom * xRandom - yRandom * yRandom)));
    }
    else if (fType == 1) {
        impactPosition = TVector3(xRandom, yRandom, -fRadius + ((xRandom * xRandom / (2 * fRadius) + yRandom * yRandom / (2 * fRadius))));
    }
    return impactPosition;
}

TVector3 TMirror::GetMirrorNormal(TVector3 impact) {
    TVector3 normal;
    if (fType == 0) {
        normal = - impact;
    }
    else if (fType == 1) {
        normal = TVector3(-impact.X(), -impact.Y(), fRadius);
    }
    return normal;
}