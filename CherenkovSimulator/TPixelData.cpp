//
//  TPixelData.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 4/6/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TPixelData.h"

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