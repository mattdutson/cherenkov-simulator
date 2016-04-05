//
//  TResponseFunction.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 4/5/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TResponseFunction_h
#define TResponseFunction_h

#include <stdio.h>
#include "TF1.h"
#include "TH1D.h"

class TResponseFunction {
    
private:
    
    Double_t fResponseTime;
    
    TF1 fResponseFucntion;
    
public:
    
    TResponseFunction(Double_t responseTime, TF1 responseFunction);
    
    TH1D ResponseHistogram(Double_t sampleTime);
    
};

#endif /* TResponseFunction_h */
