/*
 * CherenkovSimulator - TPixelData.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 *
 */

#ifndef TPixelData_h
#define TPixelData_h

#include <stdio.h>
#include "TH1D.h"

class TPixelData: public TH1D {
    
private:
    
    Double_t fX;
    
    Double_t fY;
    
public:
    
    TPixelData(Double_t x, Double_t y, TH1D histogram);
    
    void SetPosition(Double_t x, Double_t y);
    
    Double_t X();
    
    Double_t Y();
};


#endif /* TPixelData_hpp */

