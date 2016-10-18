/*
 * CherenkovSimulator - TCoordinates.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TPixelData.
 */

#include "TPixelData.hpp"

TPixelData::TPixelData(Double_t x, Double_t y, TH1D histogram): TH1D(histogram) {
    fX = x;
    fY = y;
}

void TPixelData::SetPosition(Double_t x, Double_t y) {
    fX = x;
    fY = y;
}

Double_t TPixelData::X() {
    return fX;
}

Double_t TPixelData::Y() {
    return fY;
}