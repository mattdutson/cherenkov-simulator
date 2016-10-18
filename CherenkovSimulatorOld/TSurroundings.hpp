/*
 * CherenkovSimulator - TSurroundings.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class contains information about the environment surrounding the observatory. This includes information about the ground, the atmosphere, and background noise. This class is currently being expanded.
 */

#ifndef TSurroundings_hpp
#define TSurroundings_hpp

#include "TPlane3.hpp"
#include "TScalarFunction4.hpp"
#include "TRay.hpp"
#include "TShower.hpp"

class TSurroundings {
    
private:
    
    TPlane3 fGroundPlane;
    
    TScalarFunction4* fAbsorptionCoefficients;
    
    Double_t fSpatialIntegrationStep;
    
public:
    
    TSurroundings();
    
    TSurroundings(TPlane3 groundPlane, TScalarFunction4* fAbsorptionCoefficients, Double_t spatialIntegrationStep);
    
    TPlane3 GroundPlane();
    
    Double_t GetDimmingPercentage(TShower shower, TVector3 endingPoint);
};

#endif
