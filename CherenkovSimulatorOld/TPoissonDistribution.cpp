/*
 * CherenkovSimulator - TCoordinates.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TPoissonDistribution.
 */

#include "TPoissonDistribution.hpp"

TPoissonDistribution::TPoissonDistribution(Double_t averageNumber): TF1("poisson", "[0]^x*e^(-[0])/(x!)") {
    fAverageNumber = averageNumber;
    fIntegratedDistro = std::vector<Double_t>();
    SetParameter(0, averageNumber);
    Double_t thresholdIntegral = 0.999;
    Double_t currentIntegral = 0;
    for(Int_t i = 0; currentIntegral < thresholdIntegral; i++) {
        currentIntegral += Eval(i);
        fIntegratedDistro.push_back(currentIntegral);
    }
}

Int_t TPoissonDistribution::GetNumber() {
    Double_t random = fRng->Rndm();
    Int_t i = 0;
    for (i = 0; i < fIntegratedDistro.size() && fIntegratedDistro[i] < random; i++) {}
    return i;
}