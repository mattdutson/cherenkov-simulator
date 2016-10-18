/*
 * CherenkovSimulator - TPlane3.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Represents a three-dimensional plane. This is used both in ray tracing (plane reflection) and in shower reconstruction (shower-detector plane).
 */

#ifndef TPlane3_hpp
#define TPlane3_hpp

#include "TVector3.h"

class TPlane3 {
    
private:
    
    // A vector normal to the plane
    TVector3 fNormal;
    
    // The fourth coefficient in the plane equation (ax + by + cz = d)
    Double_t fCoefficient;
    
public:
    
    /*
     * The default constructor.
     */
    TPlane3();
    
    /*
     * Creates a TPlane using the specified normal vector and some point (x, y, z) which lies on the plane.
     */
    TPlane3(TVector3 normal, TVector3 point);
    
    /*
     * Returns the plane's normal vector.
     */
    TVector3 GetNormal();
    
    /*
     * Returns the fourth coefficient in the plane equation.
     */
    Double_t GetEquationCoefficient();
    
    Double_t ShortestDistance(TVector3 point);
    
    TVector3 IntersectWithXZPlane();
    
    TVector3 ProjectOntoPlane(TVector3 point);
};

#endif
