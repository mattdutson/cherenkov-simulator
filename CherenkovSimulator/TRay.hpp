/*
 * CherenkovSimulator - TRay.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Represents a three-dimensional light ray. These rays can be propagated to objects and reflected from planes.
 */

#ifndef Ray_hpp
#define Ray_hpp

#include "TPlane3.hpp"

#include "TVector3.h"

class TRay {
    
protected:
    
    // The current time
    Double_t fTime;
    
    // The current position of the light ray
    TVector3 fPosition;
    
    // The current velocity of the light ray
    TVector3 fVelocity;
    
public:
    
    // A constant representing the speed of light
    constexpr static const Double_t fLightSpeed = 3e8;
    
    /*
     * Initializes the TRay given a position and direction.
     */
    TRay(Double_t time, TVector3 position, TVector3 direction);
    
    /*
     * Returns the current time.
     */
    Double_t GetTime();
    
    /*
     * Returns the current position of the ray.
     */
    TVector3 GetPosition();
    
    /*
     * Returns the current velocity of the ray.
     */
    TVector3 GetVelocity();
    
    /*
     * Finds the time to the specified plane.
     */
    Double_t TimeToPlane(TPlane3 plane);
    
    /*
     * Moves the ray in its current direction until it strikes a plane.
     */
    void PropagateToPlane(TPlane3 plane);
    
    /*
     * Advances the ray to the given position and advances its time using the distance between the two points.
     */
    void PropagateToPoint(TVector3 position);
    
    /*
     * Reflects the ray across a vector normal to the specified plane.
     */
    void ReflectFromPlane(TPlane3 plane);
    
    /*
     * Increments the position of the ray by making it travel in its current direction at the speed of light for the specified time.
     */
    void IncrementPosition(Double_t time);
    
    /*
     * Increments the position of the ray by making it travel in its current direction over the specified distance.
     */
    void IncrementPositionByDistance(Double_t distance);
};

#endif