/*
 * Created by Matthew Dutson on 1/18/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * A class representing a three-dimensional plane. This class is designed to be used in conjunction with the TRay and TTelescope classes.
 */

#ifndef TPlane3_h
#define TPlane3_h

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
};

#endif
