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
#include "TArrayD.h"

class TTelescope {
    
private:

    // A random object used for simulating a random ray striking the mirror
    TRandom1* fRandom = new TRandom1(12342834, 3);
    
    // 0 if the mirror is an circle, 1 if the mirror is a squre
    Int_t fMirrorType;
    
    // The radius of curvature of the mirror
    Double_t fRadius;
    
    // If the mirror is a square, this is the height. If the mirror is a circle, this is the radius.
    Double_t fSize;
    
    // The inclination angle of the telescope relative to the horizon
    Double_t fInclination;
    
    // The plane where the detection apparatus is located
    TPlane3 fFocalPlane;
    
    // The plane of the ground
    TPlane3 fGroundPlane;
    
    // A vector normal to the plane of the mirror
    TVector3 fMirrorAxis;
    
    /*
     * Simulates isotropic emission at objectPosition and detection of that radiation by the telescope.
     */
    TVector3 RayDetection(TVector3 objectPosition);
    
    /*
     * Chooses a random point on the back of the mirror based on the mirror type and dimensions.
     */
    TVector3 GetImpactPoint();
    
public:
    
    /*
     * This constructor assumes that the ground is flat and that the pixel plane is parallel to the mirror tangent plane.
     */
    TTelescope(Int_t mirrorType, Double_t radius, Double_t focalLength, Double_t inclinationAngle, Double_t size, Double_t groudHeight);
    
    /*
     * Simulates the motion of a cosmic ray shower across the field of view.
     */
    void ViewShower(TRay shower, Double_t delayTime, Int_t sampleNumber, TArrayD& yArray, TArrayD& zArray);
    
    /*
     * Simulates the detection of a single point which emits isotropically.
     */
    void ViewPoint(TVector3 position, Int_t sampleNumber, TArrayD& yArray, TArrayD& zArray);

};

#endif
