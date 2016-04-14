/*
 * CherenkovSimulator - TObservatory.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This is the overarching class which represents a cosmic ray observatory. It is capable of simulating the motion and detection of a shower, and then of reconstructing that shower. See "TObservatory.cpp" for implementation details.
 */

#ifndef TObservatory_hpp
#define TObservatory_hpp

#include "TMirror.hpp"
#include "TCoordinates.hpp"
#include "TSurroundings.hpp"
#include "TRawData.h"
#include "TSegmentedData.h"
#include "TShower.h"
#include "TCamera.hpp"
#include "TUtility.h"

class TObservatory: public TCamera {
    
private:
    
    TMirror fMirror;
    
    TCoordinates fCoordinates;
    
    TSurroundings fSurroundings;
    
    void ViewPointPrivate(TShower shower, TRawData& rawData);
    
    TPlane3 ApproximateShowerPlane(TSegmentedData data);
    
public:
    
    TObservatory(TMirror mirror, TCamera camera, TCoordinates coordinates, TSurroundings surroundings);
    
    TRawData ViewPoint(TShower shower);
    
    TRawData ViewShower(TShower shower, Double_t timeDelay);
    
    TRay ReconstructShower(TSegmentedData data);
    
};

#endif
