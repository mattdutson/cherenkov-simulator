//
//  TShower.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 1/25/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TShower.h"

TShower::TShower(TVector3 position, TVector3 direction):TRay(position, direction) {
    fRandom = new TRandom1();
}

TVector3 TShower::rayDetectionByMirror(TTelescope telescope) {
    
    Double_t heightMultiplier = fRandom->Rndm() - 0.5;
    
    Double_t widthMultiplier = fRandom->Rndm() - 0.5;
    
    Double_t z = heightMultiplier * telescope.getHeight();
    
    Double_t x = widthMultiplier * telescope.getWidth();
    
    TVector3 backPlanePoisition = *new TVector3(x, 0, z);
    
    TVector3 rotationAxis = *new TVector3(1, 0, 0);
    
    backPlanePoisition.Rotate(telescope.getInclination(), rotationAxis);
    
    TVector3 rayDirection = backPlanePoisition - fPosition;
    
    TVector3 mirrorAxis = telescope.getMirrorAxis();
    Double_t mirrorRadius = telescope.getRadius();
    
    TVector3 mirrorImpact = backPlanePoisition + mirrorAxis.Unit() * (mirrorRadius - sqrt(mirrorRadius * mirrorRadius - z * z - x * x));
    
    TRay detectedRay = *new TRay(mirrorImpact, rayDirection);
    
    detectedRay.ReflectFromPlane(*new TPlane3(-mirrorImpact, mirrorImpact));
    
    detectedRay.PropagateToPlane(telescope.getFocalPlane());
    
    return detectedRay.GetPosition();
}

