/*
 * CherenkovSimulator - TCoordinates.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains information about the telescope/observatory's orientation and offset with respect to some external frame. The methods in this class can be used to switch between reference frames.
 */

#ifndef TCoordiantes_hpp
#define TCoordiantes_hpp

#include "TVector3.h"

class TCoordinates {
    
private:
    
    Double_t fAzimuth;
    
    Double_t fElevation;
    
    TVector3 fCenterOfCurvature;
    
public:
    
    TCoordinates();
    
    TCoordinates(Double_t azimuth, Double_t elevation, TVector3 centerOfCurvature);
    
    void PositionToObservatoryFrame(TVector3& position);
    
    void PositionToExternalFrame(TVector3& position);
    
    void DirectionToObservatoryFrame(TVector3& direction);
    
    void DirectionToExternalFrame(TVector3& direction);
};

#endif
