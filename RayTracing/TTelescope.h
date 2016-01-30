/*
 * Created by Matthew Dutson on 1/29/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * A class representing a telescope which uses a spherical mirror.
 */

#ifndef TTelescope_h
#define TTelescope_h

#include "TPlane3.h"
#include "TRay.h"
#include "TGraph.h"
#include "TRandom1.h"

class TTelescope {
    
private:
    
    // The axis used for rotating the telescope by the inclination angle
    TVector3 fRotationAxis = *new TVector3(1, 0, 0);
    
    // A random object used for simulating a random ray striking the mirror
    TRandom1* fRandom = new TRandom1();
    
    // The focal length of the telescope mirror
    Double_t fFocalLength;
    
    // The radius of curvature of the mirror
    Double_t fRadius;
    
    // The height of the mirror
    Double_t fHeight;
    
    // The width of the mirror
    Double_t fWidth;
    
    // The inclination angle of the telescope relative to the horizon
    Double_t fInclination;
    
    // The plane where the detection apparatus is located
    TPlane3 fPixelPlane;
    
    // The plane of the ground
    TPlane3 fGroundPlane;
    
    // The plane tangent to the central point on the mirror
    TPlane3 fMirrorBackPlane;
    
    // A vector normal to the plane of the mirror
    TVector3 fMirrorAxis;
    
public:
    
    /*
     * This constructor assumes that the ground is flat and that the pixel plane is parallel to the mirror tangent plane.
     */
    TTelescope(Double_t focalLength, Double_t inclinationAngle, Double_t height, Double_t width, Double_t focalPlaneDistance, Double_t groudHeight);
    
    /*
     * Returns a vector normal to the mirror's tangent plane.
     */
    TVector3 GetMirrorAxis();
    
    /*
     * Returns the pixel plane of the mirror.
     */
    TPlane3 GetPixelPlane();
    
    /*
     * Returns the plane of the ground.
     */
    TPlane3 GetGroundPlane();
    
    /*
     * Simulates the motion of a cosmic ray shower across the field of view.
     */
    TGraph ViewShower(TRay shower, Double_t delayTime);

    /*
     * Simulates isotropic emission by the shower and detection of that radiation by the detector.
     */
    TVector3 RayDetectionByMirror(TRay shower);
};

#endif
