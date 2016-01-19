//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#ifndef Ray_h
#define Ray_h

#include "TVector3.h"

class TRay {
    
private:
    
    TVector3 position;
    TVector3 direction;
    
public:
    TRay(TVector3 position, TVector3 direction);
    
    void PropagateToPlane(TVector3 planeNorm, Double_t d);
    
    void ReflectFromPlane(TVector3 planeNorm);
    
    TVector3 getPosition();
};

#endif