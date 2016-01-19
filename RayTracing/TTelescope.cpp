//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#include "TRay.h"
#include "TVector3.h"

class TTelescope {
    
private:
    TVector3 fFocalPoint;
    TVector3 fCenterOfCurvature;
    TVector3 fGroundPlaneNorm;
    Double_t fGroundPlaneScalar;
    
public:
    TTelescope(TVector3 focalPoint, TVector3 centerOfCurvature, TVector3 groundPlaneNorm, Double_t groundPlaneScalar) {
        this->fFocalPoint = focalPoint;
        this->fCenterOfCurvature = centerOfCurvature;
        this->fGroundPlaneNorm = groundPlaneNorm;
        this->fGroundPlaneScalar = groundPlaneScalar;
    }
    
    TVector3 RayTrace(TVector3 pixelLocation, TVector3 mirrorImpact) {
        TRay *ray = new TRay(mirrorImpact, pixelLocation - mirrorImpact);
        TVector3 normal = fCenterOfCurvature - mirrorImpact;
        ray->ReflectFromPlane(normal);
        ray->PropagateToPlane(fGroundPlaneNorm, fGroundPlaneScalar);
        return ray->getPosition();
    }
};