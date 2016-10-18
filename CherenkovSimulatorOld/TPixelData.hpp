/*
 * CherenkovSimulator - TPixelData.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * An extension of a histogram which also includes some x and y position. This is used to associate the data from some photomultipler tube with its position.
 */

#ifndef TPixelData_hpp
#define TPixelData_hpp

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

