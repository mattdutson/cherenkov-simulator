//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#ifndef TTelescope_h
#define TTelescope_h

#include "TVector3.h"

class TTelescope {
    
private:
    TVector3 fFocalPoint;
    TVector3 fCenterOfCurvature;
    TVector3 fGroundPlaneNorm;
    Double_t fGroundPlaneScalar;
    
public:
    
    TTelescope(TVector3 focalPoint, TVector3 centerOfCurvature, TVector3 groundPlaneNorm, Double_t groundPlaneScalar);
    
    TVector3 RayTrace(TVector3 pixelLocation, TVector3 mirrorImpact);
};

#endif