//
//  TShower.hpp
//  RayTracing
//
//  Created by Matthew Dutson on 1/25/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TShower_h
#define TShower_h

#include <stdio.h>
#include "TRay.h"
#include "TVector3.h"
#include "TRandom1.h"
#include "TTelescope.h"

class TShower: TRay {
private:
    TRandom1* fRandom;
    
public:
    
    TShower(TVector3 position, TVector3 direction);
    
    TVector3 rayDetectionByMirror(TTelescope telescope);
    
};

#endif /* TShower_h */
