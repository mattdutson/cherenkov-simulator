//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#ifndef Ray_h
#define Ray_h

#include "TVector3.h"
#include "TPlane3.h"

class TRay {
    
private:
    TVector3 fPosition;
    TVector3 fDirection;
    
public:
    TRay(TVector3 position, TVector3 direction);
    void PropagateToPlane(TPlane3 plane);
    void ReflectFromPlane(TPlane3 plane);
    TVector3 GetPosition();
};

#endif