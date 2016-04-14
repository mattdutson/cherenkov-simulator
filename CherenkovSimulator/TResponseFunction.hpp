/*
 * CherenkovSimulator - TResponseFunction.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Represents the response of a photomultipler tube to a delta function. This is used when finding the output voltage for each tube.
 */

#ifndef TResponseFunction_hpp
#define TResponseFunction_hpp

#include <stdio.h>
#include "TF1.h"
#include "TH1D.h"

class TResponseFunction: public TF1 {
    
private:
    
    Double_t fResponseTime;
    
public:
    
    TResponseFunction();
    
    TResponseFunction(Double_t responseTime, TF1 responseFunction);
    
    TH1D ResponseHistogram(Double_t sampleTime);
    
    Double_t ResponseTime();
    
};

#endif /* TResponseFunction_h */
