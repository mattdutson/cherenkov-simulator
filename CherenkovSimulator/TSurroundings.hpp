/*
 * CherenkovSimulator - TSurroundings.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class contains information about the environment surrounding the observatory. This includes information about the ground, the atmosphere, and background noise. This class is currently being expanded.
 */

#ifndef TSurroundings_hpp
#define TSurroundings_hpp

#include "TPlane3.h"

class TSurroundings {
    
private:
    
    TPlane3 fGroundPlane;
    
public:
    
    TSurroundings();
    
    TSurroundings(TPlane3 groundPlane);
    
    TPlane3 GroundPlane();
};

#endif
