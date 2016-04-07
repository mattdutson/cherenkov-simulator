//
//  TUtility.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/22/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TUtility_h
#define TUtility_h

#include "TVector3.h"
#include "THistogramArray.h"

class TUtility {
public:
    
    static Double_t SumArray(std::vector<Double_t> array);
    
    static void WriteHistogramFile(TString filename, THistogramArray histograms);
};

#endif /* TUtility_h */
