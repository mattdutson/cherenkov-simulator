//
//  TTelescope.hpp
//  RayTracing
//
//  Created by Matthew Dutson on 1/25/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TTelescope_h
#define TTelescope_h

#include "TVector3.h"
#include "TPlane3.h"
#include "TGraph.h"
#include "TRandom1.h"
#include "TRay.h"

class TTelescope {
    
private:
    
    TPlane3 focalPlane;
    
    TPlane3 groundPlane;
    
    TPlane3 mirrorBackPlane;
    
    TVector3 mirrorNormal;
    
    Double_t focalLength;
    
    Double_t height;
    
    Double_t width;
    
    Double_t inclination;
    
    TRandom1* fRandom;
    
public:
    
    /*
     * The detailed constructor.
     */
    TTelescope(Double_t focalLength, Double_t inclinationAngle, Double_t height, Double_t width, Double_t focalPlaneDistance, Double_t groudHeight);
    
    Double_t getHeight();
    
    Double_t getWidth();
    
    Double_t getInclination();
    
    TVector3 getMirrorAxis();
    
    Double_t getRadius();
    
    TPlane3 getFocalPlane();
    
    TPlane3 getGroundPlane();
    
    TGraph viewShower(TRay shower, Double_t delayTime);

    TVector3 rayDetectionByMirror(TRay shower);
};

#endif /* TTelescope_h */
