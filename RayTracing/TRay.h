/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * A class representing a light ray. This ray can be reflected from or propagated to planes.
 */

#ifndef Ray_h
#define Ray_h

#include "TVector3.h"
#include "TPlane3.h"

class TRay {
    
private:
    
    // The current position of the light ray
    TVector3 fPosition;
    
    // The current direction of the light ray (this is opposite the direction the ray actually traveled)
    TVector3 fDirection;
    
public:
    
    /*
     * Initializes the TRay given a position and direction.
     */
    TRay(TVector3 position, TVector3 direction);
    
    /*
     * Moves the ray in its current direction until it strikes a plane.
     */
    void PropagateToPlane(TPlane3 plane);
    
    /*
     * Reflects the ray across a vector normal to the specified plane.
     */
    void ReflectFromPlane(TPlane3 plane);
    
    /*
     * Returns the current position of the ray.
     */
    TVector3 GetPosition();
};

#endif