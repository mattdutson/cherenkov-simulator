//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#ifndef TTelescope_h
#define TTelescope_h

#include "TVector3.h"
#include "TPlane3.h"

class TTelescope {
    
private:
    TVector3 fCenterOfCurvature;
    TPlane3 fGroundPlane;
    
public:
    TTelescope(TVector3 centerOfCurvature, TPlane3 groundPlane);
    TVector3 RayTrace(TVector3 pixelLocation, TVector3 mirrorImpact);
};

#endif