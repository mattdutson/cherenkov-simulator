/*
 * CherenkovSimulator - TUtility.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains a number of static utility methods for dealing with collections of items.
 */

#ifndef TUtility_hpp
#define TUtility_hpp

#include "THistogramList.hpp"

#include "TVector3.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"

class TUtility {
public:
    
    static Double_t SumArray(std::vector<Double_t> array);
    
    static void WriteHistogramFile(TString filename, THistogramList histograms);
    
    /*
     * Fills the histogram with the data from the input arrays.
     */
    static void FillHistogram(std::vector<Double_t> array1, std::vector<Double_t> array2, TH2D& histogram) ;
    
    /*
     * Fills the profile with the data from the input arrays.
     */
    static void FillProfile(std::vector<Double_t> array1, std::vector<Double_t> array2, TProfile& profile);
};

#endif /* TUtility_h */
