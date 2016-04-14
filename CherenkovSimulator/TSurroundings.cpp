/*
 * CherenkovSimulator - TSurroundings.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TSurroundings.
 */

#include "TSurroundings.hpp"

TSurroundings::TSurroundings() {}

TSurroundings::TSurroundings(TPlane3 groundPlane) {
    fGroundPlane = groundPlane;
}

TPlane3 TSurroundings::GroundPlane() {
    return fGroundPlane;
}