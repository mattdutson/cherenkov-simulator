/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This file contains the implementation of "TTelescope.h". See the header file for method descriptions.
 */

#include "TRay.h"
#include "TTelescope.h"

TTelescope::TTelescope(TVector3 centerOfCurvature, TPlane3 groundPlane) {
    fCenterOfCurvature = centerOfCurvature;
    fGroundPlane = groundPlane;
}

TVector3 TTelescope::RayTrace(TVector3 pixelLocation, TVector3 mirrorImpact) {
    
    // pixelLocation - mirrorImpact is a vector parallel to the ray's current trajectory
    TRay *ray = new TRay(mirrorImpact, pixelLocation - mirrorImpact);
    
    // A plane tangent to the mirror at the point where it was struck by the light ray
    TPlane3 *mirror = new TPlane3(fCenterOfCurvature - mirrorImpact, mirrorImpact);
    
    // Reflect the ray from the mirror and let it propagate to the ground
    ray->ReflectFromPlane(*mirror);
    ray->PropagateToPlane(fGroundPlane);
    return ray->GetPosition();
}