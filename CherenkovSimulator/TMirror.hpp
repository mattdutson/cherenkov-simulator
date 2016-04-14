/*
 * CherenkovSimulator - TMirror.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the optical properties of the observatory mirror. The distance where the camera is located is specified by TCamera. A TMirror instance is contained within TObservatory.
 */

#ifndef TMirror_hpp
#define TMirror_hpp

#include "TVector3.h"
#include "TRandom1.h"

class TMirror {
    
private:
    
    TRandom1* fRng = new TRandom1(12342834, 3);
    
    Int_t fType;
    
    Int_t fShape;
    
    Double_t fRadius;
    
    Double_t fSize;
    
public:
    
    TMirror();
    
    TMirror(Int_t type, Int_t shape, Double_t radius, Double_t size);
    
    Double_t Radius();
    
    TVector3 GetMirrorImpact();
    
    TVector3 GetMirrorNormal(TVector3 impact);
    
};

#endif
