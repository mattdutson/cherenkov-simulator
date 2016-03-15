//
//  TObservatory.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TObservatory_h
#define TObservatory_h

#include "TMirror.h"
#include "TCamera.h"
#include "TCoordinates.h"
#include "TSurroundings.h"
#include "TRawData.h"
#include "TSegmentedData.h"
#include "TShower.H"

class TObservatory {
    
private:
    
    TMirror fMirror;
    
    TCamera fCamera;
    
    TCoordinates fCoordinates;
    
    TSurroundings fSurroundings;
    
    void ViewPointPrivate(TShower shower, TRawData& rawData);
    
    TPlane3 ApproximateShowerPlane(TSegmentedData data);
    
public:
    
    TObservatory(TMirror mirror, TCamera camera, TCoordinates coordinates, TSurroundings surroundings);
    
    TRawData ViewPoint(TShower shower);
    
    TRawData ViewShower(TShower shower, Double_t timeDelay);
    
    TShower ReconstructShower(TSegmentedData data);
    
    TSegmentedData ParseData(TRawData data);
};

#endif /* TObservatory_h */
