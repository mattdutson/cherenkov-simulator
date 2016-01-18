//
//  Ray.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TVector3.h"
#include "TMath.h"

class Ray {
    
private:
    
    TVector3 position;
    TVector3 direction;
    
public:
    Ray(TVector3 position, TVector3 direction) {
        this->position = position;
        this->direction = direction;
    }
    
    void propagateToPlane(TVector3 planeNorm, Double_t d) {
        Double_t param = (d - planeNorm.Dot(position)) / planeNorm.Dot(direction);
        position += position * param;
    }
    
    void reflectFromPlane(TVector3 planeNorm) {
        direction.Rotate(TMath::Pi(), planeNorm);
        direction = -direction;
    }
    
    TVector3 getPosition() {
        return position;
    }
};
