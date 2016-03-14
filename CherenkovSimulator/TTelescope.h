/*
 * Created by Matthew Dutson on 1/29/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * A class representing a telescope which uses a spherical mirror.
 */

class TCamera;

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
#include "TRawData.h"
#include "TMirror.h"

class TTelescope: public TMirror {
    
private:
    
    Double_t fFocalLength;
    
    // The plane where the detection apparatus is located
    TPlane3 fFocalPlane;
    
    // The plane of the ground
    TPlane3 fGroundPlane;
    
    // The telescope camera
    TCamera* fCamera;
    
    /*
     * A private method for viewing a point which does not clear the input array.
     */
    void ViewPointPrivate(TShower shower, TRawData& data);
    
    /*
     * Simulates isotropic emission at objectPosition and detection of that radiation by the telescope.
     */
    TRay* RayDetection(TShower shower);
    
    TVector3 GetOutwardDirection(TTelescope telescope, Int_t pixel);
    
public:
    
    /*
     * This constructor makes a number of simplifying assumptions.
     */
    TTelescope(TCamera* camera, TMirror mirror, Double_t focalLength);
    
    /*
     * The detailed constructor.
     */
    TTelescope(TPlane3 groundPlane, TCamera* camera, TMirror mirror, Double_t focalLength);
    
    /*
     * Simulates the motion of a cosmic ray shower across the field of view.
     */
    void ViewShower(TShower shower, Double_t delayTime, TRawData& data);
    
    /*
     * Simulates the detection of a single point which emits isotropically.
     */
    void ViewPoint(TShower shower, TRawData& data);
    
    TVector3 GetOutwardDirection(Double_t pixelX, Double_t pixelY);
    
    /*
     * Gets the camera.
     */
    TCamera* GetCamera();
    
//    Double_t GetFocalLength();
//    
//    Double_t GetRadius();
//    
//    TVector3 GetAxis();
//
//    TVector3 GetCenterOfCurvature();
};

#endif
