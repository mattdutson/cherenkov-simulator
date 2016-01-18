//
//  Telescope.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 1/18/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

# include "Ray.cpp"

class Telescope {
    
private:
    TVector3 focalPoint;
    TVector3 centerOfCurvature;
    TVector3 groundPlaneNorm;
    Double_t groundPlaneScalar;
    
public:
    Telescope(TVector3 focalPoint, TVector3 centerOfCurvature, TVector3 groundPlaneNorm, Double_t groundPlaneScalar) {
        this->focalPoint = focalPoint;
        this->centerOfCurvature = centerOfCurvature;
        this->groundPlaneNorm = groundPlaneNorm;
        this->groundPlaneScalar = groundPlaneScalar;
    }
    
    TVector3 rayTrace(TVector3 pixelLocation, TVector3 mirrorImpact) {
        Ray ray = *new Ray(mirrorImpact, pixelLocation - mirrorImpact);
        TVector3 normal = centerOfCurvature - mirrorImpact;
        ray.reflectFromPlane(normal);
        ray.propagateToPlane(groundPlaneNorm, groundPlaneScalar);
        return ray.getPosition();
    }
};