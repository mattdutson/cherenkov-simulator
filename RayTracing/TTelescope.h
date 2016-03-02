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
#include "TMath.h"
#include "TCamera.h"
#include "TShower.h"

class TTelescope {
    
private:

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
    
    // The plane where the detection apparatus is located
    TPlane3 fFocalPlane;
    
    // The plane of the ground
    TPlane3 fGroundPlane;
    
    // A vector normal to the plane of the mirror
    TVector3 fMirrorAxis;
    
    // The telescope camera
    TCamera fCamera;
    
    /*
     * A private method for viewPoint which does not clear the input array.
     */
    void ViewPointPrivate(TShower shower, std::vector<Double_t>& xArray, std::vector<Double_t>& yArray, std::vector<Double_t>& timeArray);
    
    /*
     * Simulates isotropic emission at objectPosition and detection of that radiation by the telescope.
     */
    TRay* RayDetection(TShower shower);
    
    /*
     * Gets the mirror impact point based on the back plane impact point and the type of mirror (spherical, parabolic).
     */
    TVector3* GetMirrorImpact();
    
    /*
     * Finds a vector normal to the mirror given the point at which the ray hits the mirror.
     */
    TVector3* GetMirrorNormal(TVector3 mirrorImpact);

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
    
public:
    
    /*
     * This constructor makes a number of simplifying assumptions.
     */
    TTelescope(Short_t mirrorShape, Short_t mirrorType, Double_t radius, Double_t focalLength, Double_t fNumber, TCamera camera);
    
    /*
     * The detailed constructor.
     */
    TTelescope(Short_t mirrorShape, Short_t mirrorType, Double_t radius, Double_t focalLength, Double_t fNumber, Double_t inclination, Double_t azimuth, TVector3 centerOfCurvature, TPlane3 groundPlane, TCamera camera);
    
    /*
     * Simulates the motion of a cosmic ray shower across the field of view.
     */
    void ViewShower(TShower shower, Double_t delayTime, std::vector<Double_t>& xArray, std::vector<Double_t>& yArray, std::vector<Double_t>& timeArray);
    
    /*
     * Simulates the detection of a single point which emits isotropically.
     */
    void ViewPoint(TShower shower, std::vector<Double_t>& xArray, std::vector<Double_t>& yArray, std::vector<Double_t>& timeArray);
    
    
    /*
     * Gets the camera.
     */
    TCamera* GetCamera();

};

#endif
