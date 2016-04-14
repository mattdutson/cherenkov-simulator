/*
 * CherenkovSimulator - TAnalysis.hpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains static methods for doing statistical tests. Currently, it is used to find the spread of raw data points over the focal plane at various angles off-axis.
 */

#ifndef TAnalysis_hpp
#define TAnalysis_hpp

#include "TObservatory.hpp"
#include "TConstantIntensity.hpp"
#include "TUtility.hpp"

#include "TMath.h"
#include "TVector3.h"
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
    

};

#endif