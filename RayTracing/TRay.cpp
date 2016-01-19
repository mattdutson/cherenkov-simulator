//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.

#include "TVector3.h"
#include "TMath.h"

class TRay {
    
private:
    TVector3 fPosition;
    TVector3 fDirection;
    
public:
    TRay(TVector3 position, TVector3 direction) {
        this->fPosition = fPosition;
        this->fDirection = fDirection;
    }
    
    void PropagateToPlane(TVector3 planeNorm, Double_t d) {
        Double_t param = (d - planeNorm.Dot(fPosition)) / planeNorm.Dot(fDirection);
        fPosition += fPosition * param;
    }
    
    void ReflectFromPlane(TVector3 planeNorm) {
        fDirection.Rotate(TMath::Pi(), planeNorm);
        fDirection = -fDirection;
    }
    
    TVector3 GetPosition() {
        return fPosition;
    }
};
