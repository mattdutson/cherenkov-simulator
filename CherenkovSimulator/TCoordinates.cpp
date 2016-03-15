//
//  TCoordinates.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TCoordinates.h"

void TCoordinates::PositionToObservatoryFrame(TVector3& position) {
    position -= fCenterOfCurvature;
    DirectionToObservatoryFrame(position);
}

void TCoordinates::PositionToExternalFrame(TVector3& position) {
    DirectionToExternalFrame(position);
    position += fCenterOfCurvature;
}

void TCoordinates::DirectionToObservatoryFrame(TVector3& direction) {
    direction.RotateX(-fElevation);
    direction.RotateY(-fAzimuth);
}

void TCoordinates::DirectionToExternalFrame(TVector3& direction) {
    direction.RotateX(fElevation);
    direction.RotateY(fAzimuth);
}