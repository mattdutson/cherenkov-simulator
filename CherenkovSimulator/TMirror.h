//
//  TMirror.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TMirror_h
#define TMirror_h

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

#endif /* TMirror_hpp */
