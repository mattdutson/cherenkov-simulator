#ifndef TTelescope_h
#define TTelescope_h

#include "TVector3.h"
#include "TPlane3.h"

/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * A class representing a telescope. This telescope consists of a spherical mirror oriented with respect to some ground plane.
 */
class TTelescope {
    
private:
    
    // The position vector of the mirror's center of curvature
    TVector3 fCenterOfCurvature;
    
    // The plane of the ground
    TPlane3 fGroundPlane;
    
public:
    
    /*
     * Initializes the TTelescope given the center of curvature and the plane of the ground.
     */
    TTelescope(TVector3 centerOfCurvature, TPlane3 groundPlane);
    
    /*
     * Traces a ray detected by the telescope back to its source. pixelLocation is the position of the pixel where the ray was viewed, mirrorImpact is the location on the mirror where the ray was reflected.
     */
    TVector3 RayTrace(TVector3 pixelLocation, TVector3 mirrorImpact);
};

#endif