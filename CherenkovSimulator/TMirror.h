//
//  TMirror.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/12/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TMirror_h
#define TMirror_h

#include <stdio.h>
#include "TVector3.h"
#include "TRandom1.h"
#include "TMath.h"

#endif /* TMirror_h */

class TMirror {
    
protected:
    // A random object used for simulating a random ray striking the mirror
    TRandom1* fRandom = new TRandom1(12342834, 3);
    
    // 0 if the mirror is an circle, 1 if the mirror is a square
    Short_t fMirrorShape;
    
    // 0 if the mirror is spherical, 1 if the mirror is parabolic
    Short_t fMirrorType;
    
    // The radius of curvature of the mirror
    Double_t fRadius;
    
    // If the mirror is a square, this is the length of one side. If the mirror is a circle, this is the diameter
    Double_t fCrossDiameter;
    
    // The inclination angle of the telescope relative to the horizon
    Double_t fInclination;
    
    // The angle of the telescope relative to the positive z axis
    Double_t fAzimuth;
    
    // The location of the mirror's center of curvature
    TVector3 fCenterOfCurvature;
    
    // A vector normal to the plane of the mirror
    TVector3 fMirrorAxis;
    
public:
    
    TMirror(Short_t mirrorShape, Short_t mirrorType, Double_t radius, Double_t size, Double_t inclination, Double_t azimuth, TVector3 centerOfCurvature);
    
    TMirror(Short_t mirrorShape, Short_t mirrorType, Double_t radius, Double_t size);
    
    /*
     * Finds a vector normal to the mirror given the point at which the ray hits the mirror.
     */
    TVector3* GetMirrorNormal(TVector3 mirrorImpact);
    
    /*
     * Gets the mirror impact point based on the back plane impact point and the type of mirror (spherical, parabolic).
     */
    TVector3* GetMirrorImpact();
    
    /*
     * Rotates the vector from the telescope frame to the lab frame.
     */
    void RotateIn(TVector3& vector);
    
    /*
     * Rotates the vector from the lab frame to the telescope frame.
     */
    void RotateOut(TVector3& vector);
    
    /*
     * Translates the vector from the telescope frame to the lab frame.
     */
    void TranslateIn(TVector3& vector);
    
    /*
     * Translates the vector from the lab frame into the telescope frame.
     */
    void TranslateOut(TVector3& vector);

};
