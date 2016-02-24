/*
 * Created by Matthew Dutson on 2/23/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This file contains the implementation of "TDataPoint.h". See the header file for method descriptions.
 */

#include "TDataPoint.h"

TDataPoint::TDataPoint(Double_t time, Double_t x, Double_t y) {
    fTime = time;
    fX = x;
    fY = y;
}

Double_t TDataPoint::GetTime() {
    return fTime;
}

Double_t TDataPoint::GetX() {
    return fX;
}

Double_t TDataPoint::GetY() {
    return fY;
}