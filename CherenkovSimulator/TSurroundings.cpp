//
//  TSurroundings.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/14/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TSurroundings.h"

TSurroundings::TSurroundings() {}

TSurroundings::TSurroundings(TPlane3 groundPlane) {
    fGroundPlane = groundPlane;
}

TPlane3 TSurroundings::GroundPlane() {
    return fGroundPlane;
}