/*
 * CherenkovSimulator - TResponseFunction.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 *
 */

#ifndef TResponseFunction_h
#define TResponseFunction_h

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
