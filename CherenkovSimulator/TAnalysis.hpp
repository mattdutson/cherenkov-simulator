/*
 * CherenkovSimulator - TAnalysis.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 *
 */

#ifndef TAnalysis_h
#define TAnalysis_h

#include "TObservatory.hpp"
#include "TH2.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TArrayD.h"
#include <stdio.h>

class TAnalysis {
    
private:
    
    /*
     * Finds the RMS deviation from the average.
     */
    static Double_t FindRMSDeviation(TRawData data);
    
public:

    /*
     * Returns a graph of RMS Deviation vs angle from the mirror axis.
     */
    static void FindRMSVsAngle(std::vector<Double_t>& RMS, std::vector<Double_t>& angle, TObservatory observatory, Int_t sampleNumber, Double_t timeDelay, Double_t minAngle, Double_t maxAngle, Double_t zDistance);
    
    /*
     * Fills the histogram with the data from the input arrays.
     */
    static void FillHistogram(std::vector<Double_t> array1, std::vector<Double_t> array2, TH2D& histogram) ;
    
    /*
     * Fills the profile with the data from the input arrays.
     */
    static void FillProfile(std::vector<Double_t> array1, std::vector<Double_t> array2, TProfile& profile);
};

#endif