//
//  TOrientation.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TOrientation_h
#define TOrientation_h

#include "TVector3.h"

class TCoordinates {
    
private:
    
    Double_t fAzimuth;
    
    Double_t fElevation;
    
    TVector3 fCenterOfCurvature;
    
public:
    
    void PositionToObservatoryFrame(TVector3& position);
    
    void PositionToExternalFrame(TVector3& position);
    
    void DirectionToObservatoryFrame(TVector3& direction);
    
    void DirectionToExternalFrame(TVector3& direction);
};

#endif /* TOrientation_h */
