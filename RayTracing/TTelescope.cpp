//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#include "TRay.h"
#include "TTelescope.h"

TTelescope::TTelescope(TVector3 centerOfCurvature, TPlane3 groundPlane) {
    fCenterOfCurvature = centerOfCurvature;
    fGroundPlane = groundPlane;
}

TVector3 TTelescope::RayTrace(TVector3 pixelLocation, TVector3 mirrorImpact) {
    TRay *ray = new TRay(mirrorImpact, pixelLocation - mirrorImpact);
    TPlane3 *mirror = new TPlane3(fCenterOfCurvature - mirrorImpact, mirrorImpact);
    ray->ReflectFromPlane(*mirror);
    ray->PropagateToPlane(fGroundPlane);
    return ray->GetPosition();
}