//
//  TSurroundings.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TSurroundings_h
#define TSurroundings_h

#include "TPlane3.h"

class TSurroundings {
    
private:
    
    TPlane3 fGroundPlane;
    
public:
    
    TSurroundings();
    
    TSurroundings(TPlane3 groundPlane);
    
    TPlane3 GroundPlane();
};

#endif /* TSurroundings_h */
